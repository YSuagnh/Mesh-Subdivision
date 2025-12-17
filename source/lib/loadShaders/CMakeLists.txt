set(CMAKE_CXX_STANDARD 17)

file(GLOB lib_src 
  "*.cpp"
)

add_library(loadShaders STATIC ${lib_src})

target_include_directories(loadShaders
  PUBLIC 
  "${CMAKE_SOURCE_DIR}/include"
  "${CMAKE_SOURCE_DIR}/deps/glad/include"
)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lGL")