#include "a2.hpp"

void mesh::vertex_set_construction(const vec3* in_vertices, const vec3* in_normals, int nv, bool normals_present){
    if( vertices.size() <= nv){
        vertices.resize(nv);
    }
    for(int i = 0; i < nv; i++){
        vertices[i].world_pos = in_vertices[i];
        if(normals_present){
            vertices[i].normal = in_normals[i];
        }
    }
}