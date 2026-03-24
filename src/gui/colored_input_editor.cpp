#include "gui/colored_input_editor.h"

#include <FL/fl_draw.H>

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <unordered_set>

namespace easynmr {
namespace {

Fl_Color rgb(unsigned char r, unsigned char g, unsigned char b) {
    return fl_rgb_color(r, g, b);
}

// Style char mapping:
// A default, B carbon, C nitrogen, D oxygen, E sulfur, F phosphorus,
// G halogens, H hydrogen, I other known atoms.
Fl_Text_Display::Style_Table_Entry kStyleTable[] = {
    {rgb(74, 83, 98), FL_COURIER, 13},  // A
    {rgb(98, 106, 122), FL_COURIER, 13}, // B
    {rgb(121, 141, 220), FL_COURIER, 13}, // C
    {rgb(225, 114, 114), FL_COURIER, 13}, // D
    {rgb(222, 177, 113), FL_COURIER, 13}, // E
    {rgb(216, 154, 105), FL_COURIER, 13}, // F
    {rgb(112, 178, 143), FL_COURIER, 13}, // G
    {rgb(196, 205, 217), FL_COURIER, 13}, // H
    {rgb(151, 122, 168), FL_COURIER, 13}, // I
};

constexpr int kStyleCount = static_cast<int>(sizeof(kStyleTable) / sizeof(kStyleTable[0]));

const std::unordered_set<std::string> &known_symbols() {
    static const std::unordered_set<std::string> symbols = {
        "H",  "B",  "C",  "N",  "O",  "F",  "P",  "S",  "I",
        "Cl", "Br", "Si", "Se", "As",
        "Li", "Na", "K",  "Mg", "Ca", "Al",
        "Fe", "Co", "Ni", "Cu", "Zn", "Sn",
    };
    return symbols;
}

char style_for_aromatic(unsigned char c) {
    switch (static_cast<char>(std::tolower(c))) {
    case 'c':
        return 'B';
    case 'n':
        return 'C';
    case 'o':
        return 'D';
    case 's':
        return 'E';
    case 'p':
        return 'F';
    case 'h':
        return 'H';
    case 'b':
        return 'I';
    default:
        return 'A';
    }
}

} // namespace

ColoredInputEditor::ColoredInputEditor(int x, int y, int w, int h, const char *label)
    : Fl_Text_Editor(x, y, w, h, label),
      text_buffer_(new Fl_Text_Buffer()),
      style_buffer_(new Fl_Text_Buffer()) {
    buffer(text_buffer_);
    box(FL_DOWN_BOX);
    color(rgb(255, 255, 255));
    textfont(FL_COURIER);
    textsize(13);
    wrap_mode(Fl_Text_Display::WRAP_AT_BOUNDS, 0);
    scrollbar_width(10);

    highlight_data(style_buffer_, kStyleTable, kStyleCount, 'A', nullptr, nullptr);
    text_buffer_->add_modify_callback(on_text_modified, this);
}

ColoredInputEditor::~ColoredInputEditor() {
    if (text_buffer_ != nullptr) {
        text_buffer_->remove_modify_callback(on_text_modified, this);
    }
    delete style_buffer_;
    delete text_buffer_;
}

const char *ColoredInputEditor::value() const {
    char *raw = text_buffer_ != nullptr ? text_buffer_->text() : nullptr;
    cached_value_ = (raw != nullptr) ? raw : "";
    if (raw != nullptr) {
        std::free(raw);
    }
    return cached_value_.c_str();
}

void ColoredInputEditor::value(const char *text) {
    if (text_buffer_ == nullptr) {
        return;
    }

    suppress_change_callback_ = true;
    text_buffer_->text(text != nullptr ? text : "");
    suppress_change_callback_ = false;
    refresh_styles();
}

void ColoredInputEditor::set_syntax_mode(InputSyntaxMode mode) {
    if (syntax_mode_ == mode) {
        return;
    }
    syntax_mode_ = mode;
    refresh_styles();
}

InputSyntaxMode ColoredInputEditor::syntax_mode() const {
    return syntax_mode_;
}

void ColoredInputEditor::on_text_modified(
    int,
    int,
    int,
    int,
    const char *,
    void *cb_arg) {
    auto *self = static_cast<ColoredInputEditor *>(cb_arg);
    if (self == nullptr || self->suppress_change_callback_) {
        return;
    }
    self->refresh_styles();
    self->do_callback();
}

void ColoredInputEditor::refresh_styles() {
    if (text_buffer_ == nullptr || style_buffer_ == nullptr) {
        return;
    }

    char *raw = text_buffer_->text();
    const std::string text = (raw != nullptr) ? raw : "";
    if (raw != nullptr) {
        std::free(raw);
    }

    std::string styles = (syntax_mode_ == InputSyntaxMode::XyzLike)
        ? build_xyz_styles(text)
        : build_smiles_styles(text);

    if (styles.size() != text.size()) {
        styles.assign(text.size(), 'A');
    }
    style_buffer_->text(styles.c_str());
    redisplay_range(0, text_buffer_->length());
}

std::string ColoredInputEditor::build_smiles_styles(const std::string &text) const {
    std::string styles(text.size(), 'A');
    for (std::size_t i = 0; i < text.size(); ++i) {
        const unsigned char c = static_cast<unsigned char>(text[i]);
        if (std::isupper(c) != 0) {
            if (i + 1 < text.size()) {
                const unsigned char c2 = static_cast<unsigned char>(text[i + 1]);
                if (std::islower(c2) != 0) {
                    std::string two;
                    two.push_back(static_cast<char>(c));
                    two.push_back(static_cast<char>(c2));
                    const std::string symbol2 = normalize_symbol(two);
                    if (is_known_symbol(symbol2)) {
                        const char style = style_for_symbol(symbol2);
                        styles[i] = style;
                        styles[i + 1] = style;
                        ++i;
                        continue;
                    }
                }
            }

            const std::string symbol1(1, static_cast<char>(c));
            if (is_known_symbol(symbol1)) {
                styles[i] = style_for_symbol(symbol1);
            }
            continue;
        }

        if (std::islower(c) != 0) {
            styles[i] = style_for_aromatic(c);
        }
    }
    return styles;
}

std::string ColoredInputEditor::build_xyz_styles(const std::string &text) const {
    std::string styles(text.size(), 'A');
    std::size_t line_start = 0;
    while (line_start < text.size()) {
        std::size_t line_end = text.find('\n', line_start);
        if (line_end == std::string::npos) {
            line_end = text.size();
        }

        std::size_t token_start = line_start;
        while (token_start < line_end) {
            const unsigned char c = static_cast<unsigned char>(text[token_start]);
            if (c != ' ' && c != '\t') {
                break;
            }
            ++token_start;
        }

        std::size_t token_end = token_start;
        while (token_end < line_end && std::isalpha(static_cast<unsigned char>(text[token_end])) != 0) {
            ++token_end;
        }

        if (token_end > token_start) {
            const std::string symbol = normalize_symbol(text.substr(token_start, token_end - token_start));
            if (is_known_symbol(symbol)) {
                const char style = style_for_symbol(symbol);
                for (std::size_t i = token_start; i < token_end; ++i) {
                    styles[i] = style;
                }
            }
        }

        line_start = (line_end < text.size()) ? line_end + 1 : text.size();
    }
    return styles;
}

std::string ColoredInputEditor::normalize_symbol(const std::string &token) {
    if (token.empty()) {
        return {};
    }
    std::string symbol = token;
    symbol[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(symbol[0])));
    for (std::size_t i = 1; i < symbol.size(); ++i) {
        symbol[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(symbol[i])));
    }
    return symbol;
}

bool ColoredInputEditor::is_known_symbol(const std::string &symbol) {
    return known_symbols().count(symbol) > 0;
}

char ColoredInputEditor::style_for_symbol(const std::string &symbol) {
    if (symbol == "C") {
        return 'B';
    }
    if (symbol == "N") {
        return 'C';
    }
    if (symbol == "O") {
        return 'D';
    }
    if (symbol == "S") {
        return 'E';
    }
    if (symbol == "P") {
        return 'F';
    }
    if (symbol == "F" || symbol == "Cl" || symbol == "Br" || symbol == "I") {
        return 'G';
    }
    if (symbol == "H") {
        return 'H';
    }
    return 'I';
}

} // namespace easynmr
