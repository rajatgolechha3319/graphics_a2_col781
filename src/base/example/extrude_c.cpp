#include "../viewer/viewer.hpp"

namespace V = COL781::Viewer;
using namespace glm;

int main() {


    V::Viewer v;
    if (!v.initialize("Mesh viewer", 640, 480)) {
        return EXIT_FAILURE;
    }


    v.load_obj_file("/Users/darkelixir/Desktop/COL781/A2/my_v/meshes/cube_2.obj");
    std::vector<int> ff;
    ff.push_back(0);
    ff.push_back(1);
    ff.push_back(2);
    v.extrude_region(ff, 0.33f);
    v.view();   
}
