set(CMAKE_CXX_STANDARD 20)

file(GLOB lib_src 
  "*.cpp"
)

add_library(loadShaders STATIC ${lib_src})

target_include_directories(loadShaders
  PUBLIC 
  "${CMAKE_SOURCE_DIR}/include"
)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lGL")