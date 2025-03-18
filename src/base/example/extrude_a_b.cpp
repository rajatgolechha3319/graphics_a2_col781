#include "../viewer/viewer.hpp"

namespace V = COL781::Viewer;
using namespace glm;

int main() {

    V::Viewer v;
    if (!v.initialize("Mesh viewer", 640, 480)) {
        return EXIT_FAILURE;
    }

    // Extrude 6A
    v.create_cube_new(3,3,3);
    v.extrude(v.get_closest_face(vec3(0.0f,0.0f,0.5f)),0.33f);
    v.extrude(v.get_closest_face(vec3(0.0f,0.0f,-0.5f)),0.33f);
    v.extrude(v.get_closest_face(vec3(0.0f,0.5f,0.0f)),0.33f);
    v.extrude(v.get_closest_face(vec3(0.0f,-0.5f,0.0f)),0.33f);
    v.extrude(v.get_closest_face(vec3(0.5f,0.0f,0.0f)),0.33f);
    v.extrude(v.get_closest_face(vec3(-0.5f,0.0f,0.0f)),0.33f);


    v.view();
}
