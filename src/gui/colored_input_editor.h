#pragma once

#include <FL/Fl_Text_Buffer.H>
#include <FL/Fl_Text_Editor.H>

#include <string>

namespace easynmr {

enum class InputSyntaxMode {
    SmilesLike,
    XyzLike
};

class ColoredInputEditor : public Fl_Text_Editor {
  public:
    ColoredInputEditor(int x, int y, int w, int h, const char *label = nullptr);
    ~ColoredInputEditor() override;

    const char *value() const;
    void value(const char *text);
    void set_syntax_mode(InputSyntaxMode mode);
    InputSyntaxMode syntax_mode() const;

  private:
    static void on_text_modified(int pos, int inserted, int deleted, int restyled, const char *deleted_text, void *cb_arg);

    void refresh_styles();
    std::string build_smiles_styles(const std::string &text) const;
    std::string build_xyz_styles(const std::string &text) const;
    static std::string normalize_symbol(const std::string &token);
    static bool is_known_symbol(const std::string &symbol);
    static char style_for_symbol(const std::string &symbol);

    Fl_Text_Buffer *text_buffer_ = nullptr;
    Fl_Text_Buffer *style_buffer_ = nullptr;
    mutable std::string cached_value_;
    InputSyntaxMode syntax_mode_ = InputSyntaxMode::SmilesLike;
    bool suppress_change_callback_ = false;
};

} // namespace easynmr
