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

add_executable(catmull_clark_a src/base/example/catmull_clark_a.cpp)
target_link_libraries(catmull_clark_a viewer)

add_executable(catmull_clark_b src/base/example/catmull_clark_b.cpp)
target_link_libraries(catmull_clark_b viewer)

add_executable(catmull_clark_c src/base/example/catmull_clark_c.cpp)
target_link_libraries(catmull_clark_c viewer)

add_executable(scene src/base/example/scene.cpp)
target_link_libraries(scene viewer)

add_executable(extrude_a_b src/base/example/extrude_a_b.cpp)
target_link_libraries(extrude_a_b viewer)

add_executable(extrude_c src/base/example/extrude_c.cpp)
target_link_libraries(extrude_c viewer)

add_executable(load_obj_bunny src/base/example/load_obj_bunny.cpp)
target_link_libraries(load_obj_bunny viewer)

add_executable(load_obj_cube src/base/example/load_obj_cube.cpp)
target_link_libraries(load_obj_cube viewer)

add_executable(load_obj_spot src/base/example/load_obj_spot.cpp)
target_link_libraries(load_obj_spot viewer)

add_executable(noisy_cube src/base/example/noisy_cube.cpp)
target_link_libraries(noisy_cube viewer)

add_executable(unit_cube src/base/example/unit_cube.cpp)
target_link_libraries(unit_cube viewer)

add_executable(unit_sphere src/base/example/unit_sphere.cpp)
target_link_libraries(unit_sphere viewer)

add_executable(unit_rectangle src/base/example/unit_rectangle.cpp)
target_link_libraries(unit_rectangle viewer)

add_executable(part_8 src/a2/part_8.cpp)
target_link_libraries(part_8 glm::glm OpenGL::GL ${SDL2_LIBRARIES})
target_include_directories(part_8 PUBLIC deps/include ${SDL2_INCLUDE_DIRS})

# Output some debug info
message(STATUS "SDL2 include path: ${SDL2_INCLUDE_DIRS}")
message(STATUS "SDL2 libraries: ${SDL2_LIBRARIES}")