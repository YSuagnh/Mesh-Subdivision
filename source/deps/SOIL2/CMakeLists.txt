cmake_minimum_required(VERSION 3.4...3.8 FATAL_ERROR)

file(GLOB SOIL_src "*.c")

add_library(SOIL2 STATIC ${SOIL_src})

target_include_directories(SOIL2 PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
