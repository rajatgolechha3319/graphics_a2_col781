cmake_minimum_required(VERSION 3.0)

project(a2)
set(CMAKE_CXX_STANDARD 11)

set(GLM_INCLUDE_DIR /usr/local/Cellar/glm/1.0.1/include)
include_directories(${GLM_INCLUDE_DIR})

set(CMAKE_C_COMPILER /usr/local/bin/gcc-14)
set(CMAKE_CXX_COMPILER /usr/local/bin/g++-14)

find_package(glm REQUIRED)
find_package(OpenGL REQUIRED)
find_package(SDL2 REQUIRED)

add_library(viewer src/hw.cpp src/a2.cpp src/viewer.cpp deps/src/gl.c)
target_include_directories(viewer PUBLIC deps/include)
target_link_libraries(viewer glm::glm OpenGL::GL SDL2::SDL2)

add_executable(example src/example.cpp)
target_link_libraries(example viewer)

include_directories(/usr/local/Cellar/sdl2/2.30.11/include)