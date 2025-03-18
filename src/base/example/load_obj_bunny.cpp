#include "../viewer/viewer.hpp"

namespace V = COL781::Viewer;
using namespace glm;

int main() {

    V::Viewer v;
    if (!v.initialize("Mesh viewer", 640, 480)) {
        return EXIT_FAILURE;
    }
    v.load_obj_file("/Users/darkelixir/Desktop/COL781/A2/my_v/meshes/bunny_1k.obj");
    v.view();
}
