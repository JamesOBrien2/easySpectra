cmake_minimum_required(VERSION 3.20)
project(easySpectra VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(EASYSPECTRA_BUILD_GUI "Build FLTK desktop GUI" ON)

find_package(Python3 COMPONENTS Interpreter REQUIRED)

add_library(easynmr_core
    src/core/pipeline.cpp
    src/core/spectral_product.cpp
    src/core/spectrum.cpp
    src/core/properties.cpp
    src/core/batch_parser.cpp
)

target_include_directories(easynmr_core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_compile_definitions(easynmr_core
    PUBLIC
        EASYSPECTRA_DEFAULT_PYTHON=\"${Python3_EXECUTABLE}\"
        EASYSPECTRA_BACKEND_SCRIPT=\"${CMAKE_CURRENT_SOURCE_DIR}/backend/easynmr_backend.py\"
)

add_executable(easynmr src/cli/main.cpp)
target_link_libraries(easynmr PRIVATE easynmr_core)

add_executable(easynmr-expcheck src/cli/experimental_check.cpp)
target_link_libraries(easynmr-expcheck PRIVATE easynmr_core)

if(EASYSPECTRA_BUILD_GUI)
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
            src/gui/properties_widget.cpp
        )
        target_include_directories(easynmr-gui PRIVATE ${FLTK_INCLUDE_DIRS} ${FLTK_INCLUDE_DIR})
        if(TARGET FLTK::FLTK)
            target_link_libraries(easynmr-gui PRIVATE easynmr_core FLTK::FLTK)
        elseif(TARGET fltk::fltk)
            target_link_libraries(easynmr-gui PRIVATE easynmr_core fltk::fltk)
        else()
            target_link_libraries(easynmr-gui PRIVATE easynmr_core ${FLTK_LIBRARIES})
        endif()
        target_compile_definitions(easynmr-gui PRIVATE EASYSPECTRA_GUI_ENABLED=1)
    else()
        message(WARNING "FLTK not found: building CLI only. Set EASYSPECTRA_BUILD_GUI=OFF to silence this warning.")
    endif()
endif()
