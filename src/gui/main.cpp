#include "gui/app_window.h"

#include <FL/Fl.H>

int main() {
    Fl::lock();
    Fl::scheme("gtk+");
    Fl::background(237, 242, 248);
    Fl::background2(248, 250, 253);
    Fl::foreground(67, 77, 92);

    easynmr::AppWindow window(1320, 760, "EasyNMR Prototype");
    window.show();
    return Fl::run();
}
