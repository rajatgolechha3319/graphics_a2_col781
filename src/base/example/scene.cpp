#include "../viewer/viewer.hpp"

namespace V = COL781::Viewer;
using namespace glm;

int main() {

    V::Viewer v;
    if (!v.initialize("Mesh viewer", 640, 480)) {
        return EXIT_FAILURE;
    }

    v.load_obj_file("/Users/darkelixir/Desktop/COL781/A2/my_v/meshes/part_8.obj");
    v.catmull_clark();
    v.umbrella_update_mesh(0.01f, 100);
    v.flip_normals();
    v.view();
}
