#include <vector>
#include <set>
#include <string>
#include <map>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <glm/glm.hpp>
#define vector std::vector
#define vec3 glm::vec3

// Rules
"""
1. The edges are directed and point in one direction. The twin edge points in the opposite direction.
2. Always for a vertex choose the leftmost edge.

"""

struct vertex{
    // Add any other vertex attribs also here
    vec3 world_pos;
    vec3 normal;
    int half_edge_idx; // Should be the leftmost index, 
    // if it loops back then can be any of the half edge pointing to vertex
};

struct half_edge{
    int vertex_idx;
    int face_idx;
    int next_half_edge_idx; // Next half edge in the face in counter clockwise direction
    int twin_half_edge_idx;
};

struct face{
    int face_idx;
    int half_edge_idx; // any edge 
    vec3 face_normal;
};

class mesh{
    // Input data
    // 1. Vertices
    // 2. Vertex normals
    // 3. Faces
    // 4. Face Normals
    // 5. Boundary edges

    // Internal variables
    map<pair<int,int>, int> edge_to_half_edge; // This way we know if there is any pre-existing edge on that vertex pair that can be twin of the current edge
    vector<half_edge> half_edge_vector; // This will be used to store the half edges

    // Derived
    vector<vertex> vertices; // Derived from 1,2
    vector<half_edge> half_edges; // Derived from 3, 4, 5
    vector<face> faces; // Derived from 3, 4, 5

    // Constructors
    void vertex_set_construction(const vec3* in_vertices, const vec3* in_normals, int nv, bool normals_present);
    void face_set_construction(const vector<vector<int>> &faces, int nf);

};