#include "../viewer/viewer.hpp"

namespace V = COL781::Viewer;
using namespace glm;

int main() {

    V::Viewer v;
    if (!v.initialize("Mesh viewer", 640, 480)) {
        return EXIT_FAILURE;
    }

    v.create_cube_new(1,1,1);
    v.catmull_clark();
    v.catmull_clark();
    v.view();
    
}
