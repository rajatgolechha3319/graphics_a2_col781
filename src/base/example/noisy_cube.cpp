#include "../viewer/viewer.hpp"

namespace V = COL781::Viewer;
using namespace glm;

int main() {

    V::Viewer v;
    if (!v.initialize("Mesh viewer", 640, 480)) {
        return EXIT_FAILURE;
    }
    v.create_noisy_cube_new(10,10,10);
    // float delta = 0.01f;
    // // 10 50 100
    // int iters = 100;
    // v.umbrella_update_mesh(delta,iters);
    v.view();
}
