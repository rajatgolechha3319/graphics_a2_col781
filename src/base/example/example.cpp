#include "../viewer/viewer.hpp"

namespace V = COL781::Viewer;
using namespace glm;

int main() {

    vec3 vertices[] = {
        vec3(-0.5, -0.5, -0.5),
        vec3( 0.5, -0.5, -0.5),
        vec3(-0.5,  0.5, -0.5),
        vec3( 0.5,  0.5, -0.5),

        vec3(-0.5, -0.5,  0.5),
        vec3( 0.5, -0.5,  0.5),
        vec3(-0.5,  0.5,  0.5),
        vec3( 0.5,  0.5,  0.5),

        vec3(-0.5, -0.5, -0.5),
        vec3(-0.5, -0.5,  0.5),
        vec3(-0.5,  0.5, -0.5),
        vec3(-0.5,  0.5,  0.5),

        vec3( 0.5, -0.5, -0.5),
        vec3( 0.5, -0.5,  0.5),
        vec3( 0.5,  0.5, -0.5),
        vec3( 0.5,  0.5,  0.5),

        vec3(-0.5, -0.5, -0.5),
        vec3( 0.5, -0.5, -0.5),
        vec3(-0.5, -0.5,  0.5),
        vec3( 0.5, -0.5,  0.5),
    };

    vec3 normals[] = {
        vec3(0.0, 0.0, -1.0),
        vec3(0.0, 0.0, -1.0),
        vec3(0.0, 0.0, -1.0),
        vec3(0.0, 0.0, -1.0),

        vec3(0.0, 0.0, 1.0),
        vec3(0.0, 0.0, 1.0),
        vec3(0.0, 0.0, 1.0),
        vec3(0.0, 0.0, 1.0),

        vec3(-1.0, 0.0, 0.0),
        vec3(-1.0, 0.0, 0.0),
        vec3(-1.0, 0.0, 0.0),
        vec3(-1.0, 0.0, 0.0),

        vec3(1.0, 0.0, 0.0),
        vec3(1.0, 0.0, 0.0),
        vec3(1.0, 0.0, 0.0),
        vec3(1.0, 0.0, 0.0),

        vec3(0.0, -1.0, 0.0),
        vec3(0.0, -1.0, 0.0),
        vec3(0.0, -1.0, 0.0),
        vec3(0.0, -1.0, 0.0),
    };

    ivec3 triangles[] = {
        ivec3(0, 1, 2),
        ivec3(1, 2, 3),

        ivec3(4, 5, 6),
        ivec3(5, 6, 7),

        ivec3(8, 9, 10),
        ivec3(9, 10, 11),

        ivec3(12, 13, 14),
        ivec3(13, 14, 15),

        ivec3(16, 17, 18),
        ivec3(17, 18, 19)
    };

    ivec2 edges[] = {
        ivec2(0, 1),
        ivec2(1, 3),
        ivec2(3, 2),
        ivec2(2, 0),

        ivec2(4, 5),
        ivec2(5, 7),
        ivec2(7, 6),
        ivec2(6, 4),

        ivec2(0, 4),
        ivec2(1, 5),
        ivec2(2, 6),
        ivec2(3, 7)
    };

    V::Viewer v;
    if (!v.initialize("Mesh viewer", 640, 480)) {
        return EXIT_FAILURE;
    }

    v.load_obj_file("/Users/darkelixir/Desktop/COL781/A2/my_v/meshes/spot_control_mesh.obj");
    // v.setMesh_testing(20, 10, 12, vertices, triangles, edges, normals);
    // v.create_unit_rectangle(10, 5);
    // v.create_sphere(11, 5);
    // v.create_cube(5, 5, 5);
     v.catmull_clark();
     v.catmull_clark();
    // v.create_noisy_cube(10, 10, 10);
    // v.view(true);
    v.view();
}
