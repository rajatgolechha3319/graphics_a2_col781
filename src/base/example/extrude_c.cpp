#include "../viewer/viewer.hpp"

namespace V = COL781::Viewer;
using namespace glm;

int main() {


    V::Viewer v;
    if (!v.initialize("Mesh viewer", 640, 480)) {
        return EXIT_FAILURE;
    }


    v.create_sphere(15,10);
    std::vector<int> ff;
    ff.push_back(51);
    ff.push_back(52);
    ff.push_back(66);
    ff.push_back(67);
    v.extrude_region(ff, 0.075f);
    v.view();   
}
