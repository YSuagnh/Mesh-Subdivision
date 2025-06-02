cmake_minimum_required(VERSION 3.4...3.8 FATAL_ERROR)

add_subdirectory(glad)
add_subdirectory(glfw)
add_subdirectory(imgui)
add_subdirectory(glew-cmake)
add_subdirectory(glm)
add_subdirectory(SOIL2)

add_library(deps INTERFACE)
target_link_libraries(deps INTERFACE glad glfw imgui libglew_static glm SOIL2)
