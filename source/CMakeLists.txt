cmake_minimum_required(VERSION 3.4...3.28 FATAL_ERROR)

project(geometryDemo LANGUAGES CXX C)

option(DEBUG "load symbol table for GDB" OFF)

include(GNUInstallDirs)

add_subdirectory(deps)
add_subdirectory(lib)


# following part is a copy from glfw project
# link to glfw
# link_libraries(glfw)
include_directories("${GLFW_SOURCE_DIR}/deps")
if (MATH_LIBRARY)
  link_libraries("${MATH_LIBRARY}")
endif()

# workaround for the MS CRT deprecating parts of the standard library
if (MSVC OR CMAKE_C_SIMULATE_ID STREQUAL "MSVC")
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

if (WIN32)
  set(ICON glfw.rc)
elseif (APPLE)
  set(ICON glfw.icns)
endif()

# set(GLAD_GL "${GLFW_SOURCE_DIR}/deps/glad/gl.h")
# set(GLAD_GLES2 "${GLFW_SOURCE_DIR}/deps/glad/gles2.h")
set(GETOPT "${GLFW_SOURCE_DIR}/deps/getopt.h"
  "${GLFW_SOURCE_DIR}/deps/getopt.h")
set(TINYCTHREAD "${GLFW_SOURCE_DIR}/deps/tinycthread.h"
  "${GLFW_SOURCE_DIR}/deps/tinycthread.c")


# set compile features
add_library(main_compile_options INTERFACE)
target_compile_features(main_compile_options INTERFACE cxx_std_17)

# create main
if (DEBUG)
  add_compile_options(-g)   # GDB debug
else()
  add_compile_options(-O2)
endif()
add_executable(main main.cpp)

target_include_directories(main PRIVATE 
  "${CMAKE_SOURCE_DIR}/deps/imgui" 
)

find_package(OpenGL REQUIRED)
# find_package(GLEW REQUIRED)
target_link_libraries(main ${GLEW_LIBRARIES})

if (NOT ${OPENGL_FOUND})
  message("ERROR: OpenGL not found")
elseif(NOT ${OPENGL_GLU_FOUND})
  message("ERROR: GLU not found")
else()
  target_link_libraries(main ${OPENGL_LIBRARIES})
endif()


# target_link_libraries(loadShaders deps)
target_link_libraries(main imgui deps lib)
target_link_libraries(main main_compile_options)

# copy fonts to build dir
add_custom_target(copy_font_files ALL
  COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_SOURCE_DIR}/deps/imgui/misc/fonts 
    ${CMAKE_BINARY_DIR}/imgui/fonts
)

# copy target_link_libraries(main Threads::Threads)
if (RT_LIBRARY)
  target_link_libraries(main "${RT_LIBRARY}")
endif()

set(GUI_ONLY_BINARIES main)
# set(CONSOLE_BINARIES offscreen)

set_target_properties(${GUI_ONLY_BINARIES} ${CONSOLE_BINARIES} PROPERTIES
                      C_STANDARD 99
                      FOLDER "GLFW3/Examples")

if (MSVC)
    # Tell MSVC to use main instead of WinMain
    set_target_properties(${GUI_ONLY_BINARIES} PROPERTIES
                          LINK_FLAGS "/ENTRY:mainCRTStartup")
elseif (CMAKE_C_SIMULATE_ID STREQUAL "MSVC")
    # Tell Clang using MS CRT to use main instead of WinMain
    set_target_properties(${GUI_ONLY_BINARIES} PROPERTIES
                          LINK_FLAGS "-Wl,/entry:mainCRTStartup")
endif()

if (APPLE)
    set_target_properties(main PROPERTIES MACOSX_BUNDLE_BUNDLE_NAME "Boing")

    set_source_files_properties(glfw.icns PROPERTIES
                                MACOSX_PACKAGE_LOCATION "Resources")
    set_target_properties(${GUI_ONLY_BINARIES} PROPERTIES
                          MACOSX_BUNDLE_SHORT_VERSION_STRING ${GLFW_VERSION}
                          MACOSX_BUNDLE_LONG_VERSION_STRING ${GLFW_VERSION}
                          MACOSX_BUNDLE_ICON_FILE glfw.icns
                          MACOSX_BUNDLE_INFO_PLIST "${GLFW_SOURCE_DIR}/CMake/Info.plist.in")
endif()

