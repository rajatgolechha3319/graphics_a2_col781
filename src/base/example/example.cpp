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

//    v.load_obj_file("/Users/darkelixir/Desktop/COL781/A2/my_v/meshes/bunny_1k.obj");
//    v.load_obj_file("/Users/darkelixir/Desktop/COL781/A2/my_v/meshes/cube_2.obj");
//    std::vector<int> ff;
//    ff.push_back(0);
    // ff.push_back(1);
    // ff.push_back(2);
    // v.extrude_region(ff, 0.33f);
//    v.load_obj_file("/Users/darkelixir/Desktop/COL781/A2/my_v/meshes/spot_control_mesh.obj");
//    v.load_obj_file("/Users/darkelixir/Desktop/COL781/A2/my_v/meshes/uc.obj");

    // std::vector<int> ff;
    // ff.push_back(6);
    // ff.push_back(7);
    // ff.push_back(8);
    // ff.push_back(5);
    // // ff.push_back(2);
    // v.extrude_region(ff, 0.33f);
    // v.extrude(5, 0.3f);
    // v.extrude(6, 0.1f);
    // v.extrude(7, 0.1f);
    // v.extrude(8, 0.1f);


//    v.load_obj_file("/Users/darkelixir/Desktop/COL781/A2/my_v/meshes/dhaka.obj");
    // v.setMesh_testing(20, 10, 12, vertices, triangles, edges, normals);
    // v.create_unit_rectangle(10, 5);
    // v.create_sphere(11, 50);
    // v.create_cube_new(10,10,10);
    // v.create_noisy_cube_new(10, 10, 10);
    //  v.catmull_clark();
    // v.catmull_clark();
    // v.create_noisy_cube(10, 10, 10);
    //  v.view(true); // for umbrella


    // Extrude 6A
    // v.create_cube_new(3,3,3);
    // v.extrude(v.get_closest_face(vec3(0.0f,0.0f,0.5f)),0.33f);
    // v.extrude(v.get_closest_face(vec3(0.0f,0.0f,-0.5f)),0.33f);
    // v.extrude(v.get_closest_face(vec3(0.0f,0.5f,0.0f)),0.33f);
    // v.extrude(v.get_closest_face(vec3(0.0f,-0.5f,0.0f)),0.33f);
    // v.extrude(v.get_closest_face(vec3(0.5f,0.0f,0.0f)),0.33f);
    // v.extrude(v.get_closest_face(vec3(-0.5f,0.0f,0.0f)),0.33f);


    // v.load_obj_file("/Users/darkelixir/Desktop/COL781/A2/my_v/meshes/part_8.obj");
    // v.extrude(0,0.5f);
    v.setMesh(20, 10, 12, vertices, triangles, edges, normals);
   v.view();
}
