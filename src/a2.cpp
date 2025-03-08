#include "a2.hpp"
#define wh(i,n) while(i<n)
#define pass (void)0

// Helper functions
std::pair<int, int> get_edge(int i1, int i2) {
    return std::make_pair(std::min(i1, i2), std::max(i1, i2));
}

// Function to find edge direction in a face
// Returns true if edge is (v1->v2), false if (v2->v1)
bool get_edge_direction(const std::vector<int>& face, int v1, int v2) {
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

std::vector<int> consistent_ordering(const std::vector<std::vector<int>>& faces, int nf) {
    // This returns if faces are to be flipped or not (1 = flip, 0 = keep)
    std::vector<int> flip(nf, -1); // -1 = unvisited, 0 = keep, 1 = flip
    
    // Build edge to faces mapping
    std::map<std::pair<int, int>, std::vector<int>> edge_to_faces;
    
    // For each face, add its edges to the mapping
    int face_idx = 0;
    wh(face_idx, nf) {
        const std::vector<int>& face = faces[face_idx];
        int face_size = face.size();
        
        int i = 0;
        wh(i, face_size) {
            int v1 = face[i];
            int v2 = face[(i + 1) % face_size];
            std::pair<int, int> edge = get_edge(v1, v2);
            
            edge_to_faces[edge].push_back(face_idx);
            i++;
        }
        face_idx++;
    }
    
    // BFS traversal to propagate consistent orientation
    std::queue<int> q;
    
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
                
                const std::vector<int>& face = faces[curr_face];
                int face_size = face.size();
                
                // Process each edge of the current face
                int i = 0;
                wh(i, face_size) {
                    int v1 = face[i];
                    int v2 = face[(i + 1) % face_size];
                    std::pair<int, int> edge = get_edge(v1, v2);
                    
                    // Get actual direction of this edge in current face
                    bool curr_direction = get_edge_direction(face, v1, v2);
                    // If current face is flipped, we need to consider the reverse direction
                    if (flip[curr_face] == 1) {
                        curr_direction = !curr_direction;
                    }
                    
                    // Check all adjacent faces sharing this edge
                    std::vector<int>& adj_faces = edge_to_faces[edge];
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

std::vector<std::vector<int>> new_consistent_faces(const std::vector<std::vector<int>>& faces, int nf){
    std::vector<int> flip = consistent_ordering(faces, nf);
    std::vector<std::vector<int>> new_faces(nf);
    int face_idx = 0;
    wh(face_idx, nf) {
        const std::vector<int>& face = faces[face_idx];
        int face_size = face.size();
        
        // If face is to be flipped, reverse the order of vertices
        if (flip[face_idx] == 1) {
            for (int i = face_size - 1; i >= 0; i--) {
                new_faces[face_idx].push_back(face[i]);
            }
        } else {
            int i = 0;
            wh(i, face_size) {
                new_faces[face_idx].push_back(face[i]);
                i++;
            }
        }
        face_idx++;
    }
    return new_faces;
}

std::vector<std::vector<int>> mesh::tri_converter(const glm::ivec3* triangles, int nt){
    std::vector<std::vector<int>> faces(nt);
    int face_idx = 0;
    wh(face_idx, nt){
        faces[face_idx].push_back(triangles[face_idx].x);
        faces[face_idx].push_back(triangles[face_idx].y);
        faces[face_idx].push_back(triangles[face_idx].z);
        triangles_ivec.push_back(triangles[face_idx]); // Either tri_converter or triangulate mesh only one can be called
        face_idx++;
    }
    return faces;
}

// Mesh functions
void mesh::vertex_set_construction(const glm::vec3* in_vertices, const glm::vec3* in_normals, int nv, bool normals_present){
    if( vertices.size() <= nv){
        vertices.resize(nv);
    }
    for(int i = 0; i < nv; i++){
        vertices[i].world_pos = in_vertices[i];
        vertices_pos.push_back(in_vertices[i]);
        if(normals_present){
            vertices[i].normal = in_normals[i];
            vertices_normal.push_back(in_normals[i]);
        }
    }
}

void mesh::update_vertex (int vertex_idx, int curr_half_edge_idx){
    // This is to ensure leftmost edge is always stored in the vertex
    if(vertices[vertex_idx].half_edge_idx == -1){
        vertices[vertex_idx].half_edge_idx = curr_half_edge_idx;
    } else{
        int next_edge;
        int curr_edge = vertices[vertex_idx].half_edge_idx; // The edge point to the vertex
        if(half_edge_vector[curr_half_edge_idx].twin_half_edge_idx == -1){
            // If the new edge is not a twin of any other edge means a boundary
            // So update the vertex to point to the leftmost edge
            vertices[vertex_idx].half_edge_idx = curr_half_edge_idx;
        }
        else{
            // Iteratively go to twin and loop around the face until you get an edge pointing to vertex 
            // And do the same until either we rotate around the vertex completely then do nothing
            // Else we get an edge pointing to the vertex that is on the leftmost boundary
            int base_idx = curr_edge;
            while(true){
                next_edge = half_edge_vector[curr_edge].twin_half_edge_idx;
                if(next_edge == -1){
                    // Update the vertex to point to the leftmost edge
                    vertices[vertex_idx].half_edge_idx = curr_edge;
                    break;
                }
                // We got a twin so loop around the face to get the edge pointing to vertex
                while(half_edge_vector[next_edge].vertex_idx != vertex_idx){
                    next_edge = half_edge_vector[next_edge].next_half_edge_idx;
                }
                curr_edge = next_edge;
                if(curr_edge == base_idx){
                    // We looped around the vertex completely
                    break;
                }
            }
        }
    }
}

void mesh::triangulate_mesh(std::vector<std::vector<int>> &faces, int nf){
    // Triangulate the mesh
    std::vector<std::vector<int>> new_faces;
    int face_idx = 0;
    wh(face_idx, nf){
        const std::vector<int>& face = faces[face_idx];
        int face_size = face.size();
        if(face_size > 3){
            // Triangulate the face
            int i = 0;
            wh(i, face_size - 2){
                std::vector<int> new_face;
                new_face.push_back(face[0]);
                new_face.push_back(face[i+1]);
                new_face.push_back(face[i+2]);
                // Push to triangles ivec
                new_faces.push_back(new_face);
                triangles_ivec.push_back(glm::ivec3(face[0], face[i+1], face[i+2]));
                i++;
            }
        } else {
            new_faces.push_back(face);
            triangles_ivec.push_back(glm::ivec3(face[0], face[1], face[2]));
        }
        face_idx++;
    }
    faces = new_faces;
}
        
void mesh::face_set_construction(const std::vector<std::vector<int>> &in_faces, int nf){
    // Get the correct faces using consistent ordering
    std::vector<std::vector<int>> new_faces = new_consistent_faces(in_faces, nf);
    int face_idx = 0;
    wh(face_idx,nf){
        // Construction work
        int face_size = new_faces[face_idx].size();
        int i = 0;
        int half_edge_idx = half_edge_vector.size();
        wh(i, face_size){
            half_edge new_half_edge;
            // Edge from A to B will point B
            new_half_edge.vertex_idx = new_faces[face_idx][(i+1)%face_size];
            new_half_edge.face_idx = face_idx;
            new_half_edge.next_half_edge_idx = half_edge_idx + (i+1)%face_size;

            std::pair<int, int> twin_edge = std::make_pair(new_faces[face_idx][(i+1)%face_size], new_faces[face_idx][i]);
            std::pair<int, int> curr_edge = std::make_pair(new_faces[face_idx][i], new_faces[face_idx][(i+1)%face_size]);

            // Check if twin is in the map
            if(edge_to_half_edge.find(twin_edge) != edge_to_half_edge.end()){
                int twin_idx = edge_to_half_edge[twin_edge];
                new_half_edge.twin_half_edge_idx = twin_idx;
                half_edge_vector[twin_idx].twin_half_edge_idx = half_edge_idx + i;
            } else {
                new_half_edge.twin_half_edge_idx = -1;
                edge_to_half_edge[curr_edge] = half_edge_idx + i;
            }
            half_edge_vector.push_back(new_half_edge);
            update_vertex(new_faces[face_idx][i], half_edge_idx + i);
            boundary_edges.insert(get_edge(new_faces[face_idx][i], new_faces[face_idx][(i+1)%face_size]));
            i++;
        }
        // Update faces
        face new_face;
        new_face.face_idx = face_idx;
        new_face.half_edge_idx = half_edge_idx;
        // Skip normal
        faces.push_back(new_face);
        face_idx++;
    }
    // Convert boundary edges to vector
    for(const auto& edge_ : boundary_edges){
        boundary_edges_vec.emplace_back(edge_.first, edge_.second);
    }
}

void mesh::print_ds() {
    std::cout << "===== Mesh Data Structure =====" << std::endl;

    // Print vertices
    std::cout << "\nVertices:" << std::endl;
    for (size_t i = 0; i < vertices.size(); ++i) {
        const vertex& v = vertices[i];
        std::cout << "Vertex " << i << ": "
                  << "World Pos = (" << v.world_pos.x << ", " << v.world_pos.y << ", " << v.world_pos.z << "), "
                  << "Normal = (" << v.normal.x << ", " << v.normal.y << ", " << v.normal.z << "), "
                  << "Half-Edge Index = " << v.half_edge_idx << std::endl;
    }

    // Print half-edges
    std::cout << "\nHalf-Edges:" << std::endl;
    for (size_t i = 0; i < half_edge_vector.size(); ++i) {
        const half_edge& he = half_edge_vector[i];
        std::cout << "Half-Edge " << i << ": "
                  << "Vertex Index = " << he.vertex_idx << ", "
                  << "Face Index = " << he.face_idx << ", "
                  << "Next Half-Edge Index = " << he.next_half_edge_idx << ", "
                  << "Twin Half-Edge Index = " << he.twin_half_edge_idx << std::endl;
    }

    // Print faces
    std::cout << "\nFaces:" << std::endl;
    for (size_t i = 0; i < faces.size(); ++i) {
        const face& f = faces[i];
        std::cout << "Face " << i << ": "
                  << "Half-Edge Index = " << f.half_edge_idx << ", "
                  << "Face Normal = (" << f.face_normal.x << ", " << f.face_normal.y << ", " << f.face_normal.z << ")" << std::endl;
    }

    // Print edge-to-half-edge mapping
    std::cout << "\nEdge-to-Half-Edge Mapping:" << std::endl;
    for (const auto& entry : edge_to_half_edge) {
        const std::pair<int, int>& edge = entry.first;
        int half_edge_idx = entry.second;
        std::cout << "Edge (" << edge.first << ", " << edge.second << ") -> Half-Edge Index = " << half_edge_idx << std::endl;
    }

    std::cout << "===== End of Mesh Data Structure =====" << std::endl;
}

// Getter Function Implementations
glm::vec3* mesh::get_vertices_pos(){
    return vertices_pos.data();
}

glm::vec3* mesh::get_vertices_normal(){
    return vertices_normal.data();
}

glm::ivec3* mesh::get_triangles(){
    return triangles_ivec.data();
}

glm::ivec2* mesh::get_boundary_edges(){
    return boundary_edges_vec.data();
}

int mesh::get_num_vertices(){
    return vertices.size();
}

int mesh::get_num_faces(){
    return faces.size();
}

int mesh::get_num_boundary_edges(){
    return boundary_edges_vec.size();
}
