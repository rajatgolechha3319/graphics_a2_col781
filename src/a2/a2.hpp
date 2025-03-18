#include <vector>
#include <set>
#include <string>
#include <map>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <fstream>
#include <glm/glm.hpp>

#ifndef A2_HPP
#define A2_HPP
// Rules

// 1. The edges are directed and point in one direction. The twin edge points in the opposite direction.
// 2. Always for a vertex choose the leftmost edge.



struct vertex{
    // Add any other vertex attribs also here
    glm::vec3 world_pos;
    glm::vec3 normal;
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
    glm::vec3 face_normal;
    float area;
};

class mesh{
    public:
    // Input data
    // 1. Vertices
    // 2. Vertex normals
    // 3. Faces
    // 4. Face Normals
    // 5. Boundary edges

    // Internal variables
    std::map<std::pair<int,int>, int> edge_to_half_edge; // This way we know if there is any pre-existing edge on that vertex pair that can be twin of the current edge
    std::set<std::pair<int,int>> boundary_edges; // This will be used to store the boundary edges
    std::map<std::pair<int,int>, std::pair<int,int>> boundary_edge_to_face; // This will be used to store the faces corresponding to the boundary edges
    std::vector<glm::ivec2> boundary_edges_vec; // This will be used to store the boundary edges in vector form
    std::vector<std::vector<int>> faces_store; // This will be used to store the faces in vector form before triangulation
    std::map<int,int> face_idx_tri_to_poly; // This will be used to store the mapping from triangle face index to polygon face index
    std::map<int,int> face_idx_poly_to_tri; // Call constructor map_maker to make this
    
    // For getters
    std::vector<glm::vec3> vertices_pos; // Derived
    std::vector<glm::vec3> vertices_normal; // Derived
    std::vector<glm::ivec3> triangles_ivec; // Derived
    // Also boundary_edges_vec declared above
    
    
    // Derived
    std::vector<vertex> vertices; // Derived from 1,2
    std::vector<half_edge> half_edge_vector; // This will be used to store the half edges
    // std::vector<half_edge> half_edges; // Derived from 3, 4, 5 // We are using half_edge_vector for this
    std::vector<face> faces; // Derived from 3, 4, 5

    // Helper functions
    void update_vertex(int vertex_idx, int curr_half_edge_idx);
    void face_normal_gen(int idx);
    void vertex_normal_update(int idx);
    glm::vec3 get_umbrella_del(int idx);
    void umbrella_update_all(float delta, int iters);
    std::pair<std::vector<glm::vec3>, std::vector<std::vector<int>>>  catmull_clark_subdivision();
    void vertex_normal_update_mode_all();
    std::pair<std::vector<glm::vec3>, std::vector<std::vector<int>>> extrude_face(int face_idx, float d);
    std::pair<std::vector<glm::vec3>, std::vector<std::vector<int>>> extrude_region(std::vector<int> faces_, float d);
    float get_dist(glm::vec3 a, glm::vec3 b, glm::vec3 c, glm::vec3 p);
    int get_closest_face(glm::vec3 p);


    // Constructors
    std::vector<std::vector<int>> tri_converter(const glm::ivec3* triangles, int nt);
    void vertex_set_construction(const glm::vec3* in_vertices, const glm::vec3* in_normals, int nv, bool normals_present);
    void triangulate_mesh(std::vector<std::vector<int>> &faces, int nf);
    void face_set_construction(const std::vector<std::vector<int>> &faces, int nf);
    void map_maker();

    // Debugging functions
    void print_ds();
    void clear_mesh();

    // Getters
    glm::vec3* get_vertices_pos();
    glm::vec3* get_vertices_normal();
    glm::ivec3* get_triangles();
    glm::ivec2* get_boundary_edges();
    int get_num_vertices();
    int get_num_faces();
    int get_num_boundary_edges();

};

#endif