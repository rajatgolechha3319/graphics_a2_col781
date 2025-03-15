cmake_minimum_required(VERSION 3.10)
project(a2)
set(CMAKE_CXX_STANDARD 11)

# Set architecture to match your SDL2 library
set(CMAKE_OSX_ARCHITECTURES "x86_64")

# Set compilers
set(CMAKE_C_COMPILER /usr/local/bin/gcc-14)
set(CMAKE_CXX_COMPILER /usr/local/bin/g++-14)

# Find packages
find_package(OpenGL REQUIRED)

# Manually specify SDL2 paths since there might be architecture issues
set(SDL2_PATH "/usr/local/Cellar/sdl2/2.30.11")
set(SDL2_INCLUDE_DIRS "${SDL2_PATH}/include")
set(SDL2_LIBRARIES "${SDL2_PATH}/lib/libSDL2.dylib")

# GLM setup
set(GLM_INCLUDE_DIR "/usr/local/Cellar/glm/1.0.1/include")
include_directories(${GLM_INCLUDE_DIR})
find_package(glm REQUIRED)

# Include directories
include_directories(
    ${SDL2_INCLUDE_DIRS}
    deps/include
)

# Create library
add_library(viewer src/base/hw/hw.cpp src/a2/a2.cpp src/base/viewer/viewer.cpp deps/src/gl.c)
target_include_directories(viewer PUBLIC deps/include ${SDL2_INCLUDE_DIRS})
target_link_libraries(viewer glm::glm OpenGL::GL ${SDL2_LIBRARIES})

# Create executable
add_executable(example src/base/example/example.cpp)
target_link_libraries(example viewer)

# Output some debug info
message(STATUS "SDL2 include path: ${SDL2_INCLUDE_DIRS}")
message(STATUS "SDL2 libraries: ${SDL2_LIBRARIES}")