set(CMAKE_CXX_STANDARD 11)

file(GLOB imgui_src
  "*.cpp"
)

add_library(imgui ${imgui_src})

target_include_directories(imgui PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/backends)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lGL")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lGL")

target_link_libraries(imgui PRIVATE glfw)
