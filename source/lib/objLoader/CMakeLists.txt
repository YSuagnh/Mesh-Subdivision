cmake_minimum_required(VERSION 3.4...3.28 FATAL_ERROR)

add_library(objLoader STATIC "trackball.cc")

target_include_directories(objLoader INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}"
)