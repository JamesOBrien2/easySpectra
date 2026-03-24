cmake_minimum_required(VERSION 3.20)
project(EasyNMR VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(EASYNMR_BUILD_GUI "Build FLTK desktop GUI" ON)

find_package(Python3 COMPONENTS Interpreter REQUIRED)

add_library(easynmr_core
    src/core/pipeline.cpp
    src/core/spectral_product.cpp
    src/core/spectrum.cpp
)

target_include_directories(easynmr_core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_compile_definitions(easynmr_core
    PUBLIC
        EASYNMR_DEFAULT_PYTHON=\"${Python3_EXECUTABLE}\"
        EASYNMR_BACKEND_SCRIPT=\"${CMAKE_CURRENT_SOURCE_DIR}/backend/easynmr_backend.py\"
)

add_executable(easynmr src/cli/main.cpp)
target_link_libraries(easynmr PRIVATE easynmr_core)

if(EASYNMR_BUILD_GUI)
    find_package(FLTK QUIET)
    if(FLTK_FOUND)
        if(DEFINED FLTK_USE_FILE AND EXISTS "${FLTK_USE_FILE}")
            include("${FLTK_USE_FILE}")
        endif()
        add_executable(easynmr-gui
            src/gui/main.cpp
            src/gui/app_window.cpp
            src/gui/colored_input_editor.cpp
            src/gui/spectrum_widget.cpp
            src/gui/structure_widget.cpp
            src/gui/workflow_progress_widget.cpp
        )
        target_include_directories(easynmr-gui PRIVATE ${FLTK_INCLUDE_DIRS} ${FLTK_INCLUDE_DIR})
        if(TARGET FLTK::FLTK)
            target_link_libraries(easynmr-gui PRIVATE easynmr_core FLTK::FLTK)
        elseif(TARGET fltk::fltk)
            target_link_libraries(easynmr-gui PRIVATE easynmr_core fltk::fltk)
        else()
            target_link_libraries(easynmr-gui PRIVATE easynmr_core ${FLTK_LIBRARIES})
        endif()
        target_compile_definitions(easynmr-gui PRIVATE EASYNMR_GUI_ENABLED=1)
    else()
        message(WARNING "FLTK not found: building CLI only. Set EASYNMR_BUILD_GUI=OFF to silence this warning.")
    endif()
endif()
