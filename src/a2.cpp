#include "a2.hpp"
#define wh(i,n) while(i<n)

// Helper functions
pair<int, int> get_edge(int i1, int i2) {
    return make_pair(min(i1, i2), max(i1, i2));
}

// Function to find edge direction in a face
// Returns true if edge is (v1->v2), false if (v2->v1)
bool get_edge_direction(const vector<int>& face, int v1, int v2) {
    int size = face.size();
    int i = 0;
    wh(i, size) {
        if (face[i] == v1 && face[(i + 1) % size] == v2) {
            return true;
        }
        i++;
    }
    return false;
}

vector<int> consistent_ordering(const vector<vector<int>>& faces, int nf) {
    // This returns if faces are to be flipped or not (1 = flip, 0 = keep)
    vector<int> flip(nf, -1); // -1 = unvisited, 0 = keep, 1 = flip
    
    // Build edge to faces mapping
    map<pair<int, int>, vector<int>> edge_to_faces;
    
    // For each face, add its edges to the mapping
    int face_idx = 0;
    wh(face_idx, nf) {
        const vector<int>& face = faces[face_idx];
        int face_size = face.size();
        
        int i = 0;
        wh(i, face_size) {
            int v1 = face[i];
            int v2 = face[(i + 1) % face_size];
            pair<int, int> edge = get_edge(v1, v2);
            
            edge_to_faces[edge].push_back(face_idx);
            i++;
        }
        face_idx++;
    }
    
    // BFS traversal to propagate consistent orientation
    queue<int> q;
    
    // Process each connected component
    int start_face = 0;
    wh(start_face, nf) {
        if (flip[start_face] == -1) {
            // Start a new component with this face
            q.push(start_face);
            flip[start_face] = 0; // Choose not to flip the first face in each component
            
            while (!q.empty()) {
                int curr_face = q.front();
                q.pop();
                
                const vector<int>& face = faces[curr_face];
                int face_size = face.size();
                
                // Process each edge of the current face
                int i = 0;
                wh(i, face_size) {
                    int v1 = face[i];
                    int v2 = face[(i + 1) % face_size];
                    pair<int, int> edge = get_edge(v1, v2);
                    
                    // Get actual direction of this edge in current face
                    bool curr_direction = get_edge_direction(face, v1, v2);
                    // If current face is flipped, we need to consider the reverse direction
                    if (flip[curr_face] == 1) {
                        curr_direction = !curr_direction;
                    }
                    
                    // Check all adjacent faces sharing this edge
                    vector<int>& adj_faces = edge_to_faces[edge];
                    int adj_idx = 0;
                    wh(adj_idx, adj_faces.size()) {
                        int adj_face = adj_faces[adj_idx];
                        if (adj_face != curr_face && flip[adj_face] == -1) {
                            // Get direction of the edge in adjacent face
                            bool adj_direction = get_edge_direction(faces[adj_face], v1, v2);
                            
                            // If directions match, flip the adjacent face
                            // If directions are opposite, keep the adjacent face as is
                            if (curr_direction == adj_direction) {
                                flip[adj_face] = 1; // Flip
                            } else {
                                flip[adj_face] = 0; // Keep
                            }
                            
                            // Add to queue for further processing
                            q.push(adj_face);
                        }
                        adj_idx++;
                    }
                    i++;
                }
            }
        }
        start_face++;
    }
    
    return flip;
}


// Mesh functions
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

void mesh::face_set_construction(const vector<vector<int>> &faces, int nf){
    // For each face we have a set of edges

}