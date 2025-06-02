cmake_minimum_required(VERSION 3.4...3.8 FATAL_ERROR)

add_subdirectory(objLoader)
add_subdirectory(loadShaders)
add_subdirectory(algorithms)

add_library(lib INTERFACE)
target_link_libraries(lib INTERFACE objLoader loadShaders subdiv)
