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
        // Boundary Edge dealing
        boundary_edges.insert(std::make_pair(triangles[face_idx].x, triangles[face_idx].y));
        boundary_edges.insert(std::make_pair(triangles[face_idx].y, triangles[face_idx].z));
        boundary_edges.insert(std::make_pair(triangles[face_idx].z, triangles[face_idx].x));
        face_idx++;
    }
    for(const auto& edge_ : boundary_edges){
        boundary_edges_vec.emplace_back(edge_.first, edge_.second);
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
    // std::cout << "Updating vertex " << vertex_idx << std::endl;
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
    // std::cout << "Updated vertex " << vertex_idx << std::endl;


    // Print the half-edge vector
    // for(int i = 0; i < half_edge_vector.size(); i++){
    //     std::cout << "Half-Edge " << i << ": " << half_edge_vector[i].vertex_idx << " " << half_edge_vector[i].face_idx << " " << half_edge_vector[i].next_half_edge_idx << " " << half_edge_vector[i].twin_half_edge_idx << std::endl;
    // }

}

void mesh::face_normal_gen(int idx){
    // First fetch the face
    face curr_face = faces[idx];
    // Now generate normal using the direction of vector area
    glm::vec3 res = glm::vec3(0.0f,0.0f,0.0f);
    // Get v0
    glm::vec3 v0 = vertices[half_edge_vector[curr_face.half_edge_idx].vertex_idx].world_pos;
    int base_id = half_edge_vector[curr_face.half_edge_idx].vertex_idx;
    half_edge curr_edge = half_edge_vector[half_edge_vector[curr_face.half_edge_idx].next_half_edge_idx];
    while(curr_edge.vertex_idx != base_id){
        res += glm::cross(vertices[curr_edge.vertex_idx].world_pos - v0, vertices[half_edge_vector[curr_edge.next_half_edge_idx].vertex_idx].world_pos - v0);
        // Current = next
        curr_edge = half_edge_vector[curr_edge.next_half_edge_idx];
    }
    faces[idx].area = std::abs(glm::length(res) * 0.5f);
    res = glm::normalize(res);
    // Update in face
    faces[idx].face_normal = res;

}

void mesh::vertex_normal_update(int idx){
    // If vertex normal is not given it is set to 0,0,0 then run update else return
    float eps = 1.42e-5; // Random ass number
    if(glm::length(vertices[idx].normal) > eps){
        return;
    }

    // // We will use area weighted normal approach
    // float area_sum = 0.0f;
    // glm::vec3 weighted_normal = glm::vec3(0.0f,0.0f,0.0f);
    // // Leftmost half edge of vertex
    // half_edge v_edge = half_edge_vector[vertices[idx].half_edge_idx];
    // int base_face = v_edge.face_idx;
    // // Add details of this face
    // area_sum += faces[v_edge.face_idx].area;
    // weighted_normal += faces[v_edge.face_idx].area * faces[v_edge.face_idx].face_normal;
    // v_edge = half_edge_vector[v_edge.next_half_edge_idx];
    // if(v_edge.twin_half_edge_idx != -1){
    //     v_edge = half_edge_vector[v_edge.twin_half_edge_idx];
    //     while(v_edge.face_idx != base_face){
    //         // Add details
    //         area_sum += faces[v_edge.face_idx].area;
    //         weighted_normal += faces[v_edge.face_idx].area * faces[v_edge.face_idx].face_normal;
    //         v_edge = half_edge_vector[v_edge.next_half_edge_idx];
    //         if(v_edge.twin_half_edge_idx != -1){break;}
    //         v_edge = half_edge_vector[v_edge.twin_half_edge_idx];
    //     }
    //     weighted_normal = weighted_normal / area_sum;
    //     vertices[idx].normal = glm::normalize(weighted_normal);
    // }
    // else{
    //     // Update vertex normal to this normal
    //     vertices[idx].normal = glm::normalize(weighted_normal);
    // }

    // std::cout << " Called for vertex " << idx << '\n';
    float area_sum = 0.0f;
    glm::vec3 weighted_normal(0.0f,0.0f,0.0f);
    int he_idx = vertices[idx].half_edge_idx;
    int curr_idx = he_idx;
    int face_idx;
    half_edge current_edge;
    do{
        // Updates for current face
        current_edge = half_edge_vector[curr_idx];
        face_idx = current_edge.face_idx;

        area_sum += faces[face_idx].area;
        weighted_normal += faces[face_idx].area * faces[face_idx].face_normal;

        // Go to next
        current_edge = half_edge_vector[current_edge.next_half_edge_idx];
        // Go to twin if exist
        if(current_edge.twin_half_edge_idx == -1){break;}
        curr_idx = current_edge.twin_half_edge_idx;
        current_edge = half_edge_vector[current_edge.twin_half_edge_idx];

    } while(curr_idx != he_idx);

    weighted_normal /= area_sum;
    vertices[idx].normal = glm::normalize(weighted_normal);
    vertices_normal[idx] = vertices[idx].normal;
    

}

void mesh::vertex_normal_update_mode_all(){
    int v_iter = 0;
    int v_lim = vertices.size();
    wh(v_iter, v_lim){
        vertex_normal_update(v_iter);
        v_iter++;
    }
}

glm::vec3 mesh::get_umbrella_del(int idx) {
    half_edge curr_edge = half_edge_vector[vertices[idx].half_edge_idx];
    glm::vec3 res = glm::vec3(0.0f);
    int he_idx = vertices[idx].half_edge_idx;
    int curr_idx = he_idx;
    int count = 0;

    do {
        // Move to the next edge in the loop
        curr_idx = curr_edge.next_half_edge_idx;
        curr_edge = half_edge_vector[curr_idx];

        // Add neighbor contribution
        res += vertices[curr_edge.vertex_idx].world_pos - vertices[idx].world_pos;
        count++;

        // Handle boundary case
        if (curr_edge.twin_half_edge_idx == -1) {
            // Traverse back to find the last boundary neighbor
            int last_idx = he_idx;
            half_edge last_edge = half_edge_vector[last_idx];
            while (half_edge_vector[last_edge.next_half_edge_idx].vertex_idx != idx) {
                last_idx = last_edge.next_half_edge_idx;
                last_edge = half_edge_vector[last_idx];
            }
            res += vertices[last_edge.vertex_idx].world_pos - vertices[idx].world_pos;
            count++;
            break;
        }

        // Move to the twin edge
        curr_idx = curr_edge.twin_half_edge_idx;
        curr_edge = half_edge_vector[curr_idx];

    } while (curr_idx != he_idx);

    return count > 0 ? res / (float)count : glm::vec3(0.0f);
}

void mesh::umbrella_update_all(float delta, int iters) {

    // Updated
    int iter_loop = 0;
    int v_iter = 0;
    int v_lim = vertices.size();
    while(iter_loop < iters){
        v_iter = 0;
        wh(v_iter, v_lim) {
            vertices[v_iter].world_pos += delta * get_umbrella_del(v_iter);
            vertices_pos[v_iter] = vertices[v_iter].world_pos;
            v_iter++;
        }
        iter_loop++;
    }
    // Update vertex normals
    v_iter = 0;
    wh(v_iter, v_lim) {
        vertex_normal_update(v_iter);
        v_iter++;
    }
    
}

std::pair<std::vector<glm::vec3>, std::vector<std::vector<int>>> mesh::catmull_clark_subdivision(){
    // First we will update the vertices of existing mesh
    // Step 1 : Compute Face Points for each face
    std::vector<glm::vec3> face_points;
    // We will be using faces store for this

    // Construce boundary_edge_to_face
    std::map<std::pair<int,int>,int> vis;
    int iter1 = 0;
    int iter2 = 0;
    wh(iter1,faces_store.size()){
        auto face = faces_store[iter1];
        glm::vec3 face_point = glm::vec3(0.0f);
        for(auto vertex : face){
            face_point += vertices[vertex].world_pos;
        }
        // Iterate over all the edges and add this face as the first if not visited else second
        int n_f = face.size();
        iter2 = 0;
        wh(iter2,n_f){
            auto edge_key = get_edge(face[iter2], face[(iter2+1)%n_f]);
            if (vis.count(edge_key)) {
                boundary_edge_to_face[edge_key].second = iter1;
            } else {
                boundary_edge_to_face[edge_key] = std::make_pair(iter1, -1);
                vis[edge_key] = 42;
            }
            iter2++;
        }
        face_point /= (float)face.size();
        face_points.push_back(face_point);
        iter1++;
    }

    // Print the face points with index
    iter1 = 0;
    wh(iter1, face_points.size()){
        std::cout << iter1 << " " << face_points[iter1].x << " " << face_points[iter1].y << " " << face_points[iter1].z << "\n";
        iter1++;
    }


    // Step 2 : Compute Edge Points for each edge
    // Need to do this only for boundary edges
    std::map<std::pair<int,int>, glm::vec3> edge_points;
    for(auto edge : boundary_edges){
        glm::vec3 edge_point = 0.25f * (
            vertices[edge.first].world_pos +
            vertices[edge.second].world_pos +
            // Face 1
            face_points[boundary_edge_to_face[get_edge(edge.first, edge.second)].first] +
            // Face 2
            face_points[boundary_edge_to_face[get_edge(edge.first, edge.second)].second]
        );
        // std::cout << "Computation for edge " << edge.first << " " << edge.second << " \n";
        // std::cout << " V1 " << vertices[edge.first].world_pos.x << " " << vertices[edge.first].world_pos.y << " " << vertices[edge.first].world_pos.z << "\n";
        // std::cout << " V2 " << vertices[edge.second].world_pos.x << " " << vertices[edge.second].world_pos.y << " " << vertices[edge.second].world_pos.z << "\n";
        // std::cout << " F1 " << face_points[boundary_edge_to_face[get_edge(edge.first, edge.second)].first].x << " " << face_points[boundary_edge_to_face[get_edge(edge.first, edge.second)].first].y << " " << face_points[boundary_edge_to_face[get_edge(edge.first, edge.second)].first].z << "\n";
        // std::cout << " F2 " << face_points[boundary_edge_to_face[get_edge(edge.first, edge.second)].second].x << " " << face_points[boundary_edge_to_face[get_edge(edge.first, edge.second)].second].y << " " << face_points[boundary_edge_to_face[get_edge(edge.first, edge.second)].second].z << "\n";
        edge_points[get_edge(edge.first,edge.second)] = edge_point;
    }

    for(auto edge_ : edge_points){
        std::cout << edge_.first.first << " " << edge_.first.second << " " << edge_.second.x << " " << edge_.second.y << " " << edge_.second.z << "\n";
    }

    // Step 3 : Compute Vertex Points for each vertex
    // We will do this in a fast way iterating over the faces and creating a sum total for each vertex
    int v = get_num_vertices();
    std::vector<glm::vec3> vertice_face_sum(v, glm::vec3(0.0f));
    std::vector<int> vertice_face_count(v, 0);
    std::vector<glm::vec3> vertice_edge_sum(v, glm::vec3(0.0f));
    std::vector<int> vertice_edge_count(v, 0);
    int face_size;
    int n_faces = faces_store.size();
    int j = 0;
    int i = 0;
    wh(j, n_faces){
        // These are the original faces
        face_size = faces_store[j].size();
        i = 0;
        wh(i, face_size){
            // Add face point
            vertice_face_sum[faces_store[j][i]] += face_points[j];
            vertice_face_count[faces_store[j][i]]++;
            i++;
        }
        j++;
    }

    // Now go over edge points map to update the edge points
    for(auto edge : edge_points){
        vertice_edge_sum[edge.first.first] += edge.second;
        vertice_edge_sum[edge.first.second] += edge.second;
        vertice_edge_count[edge.first.first]++;
        vertice_edge_count[edge.first.second]++;
    }

    // Now update the vertices
    std::vector<glm::vec3> new_vertices;
    j = 0;
    wh(j, v){
        // Compute the new vertex
        glm::vec3 new_vertex = (
            ((vertice_face_sum[j]) / (float)(vertice_face_count[j]) +
            (vertice_edge_sum[j]) / (float)(vertice_edge_count[j]) * 2.0f +
            (vertices[j].world_pos) * (float)(vertice_face_count[j] - 3.0f)) / (float)(vertice_face_count[j])
        );
        new_vertices.push_back(new_vertex);
        j++;
    }
    // Add the edge points to the new vertices
    std::map<std::pair<int,int>, int> edge_point_idx;
    for(auto edge : edge_points){
        edge_point_idx[edge.first] = new_vertices.size();
        new_vertices.push_back(edge.second);
    }
    // Add the face points to the new vertices
    std::map<int, int> face_point_idx;
    j = 0;
    wh(j, n_faces){
        face_point_idx[j] = new_vertices.size();
        new_vertices.push_back(face_points[j]);
        j++;
    }

    std::vector<std::vector<int>> new_faces;
    // Process each face for the new mesh
    j = 0;
    wh(j, n_faces){
        // These were the original faces
        face_size = faces_store[j].size();
        i = 0;
        wh(i, face_size){
            std::vector<int> new_face;
            new_face.push_back(faces_store[j][i]);
            // Edge point from [j][i] to [j][(i+1)%face_size]
            new_face.push_back(edge_point_idx[get_edge(faces_store[j][i], faces_store[j][(i+1)%face_size])]);
            // Face point of face j
            new_face.push_back(face_point_idx[j]);
            // Edge point from [j][(i+face_size-1)%face_size] to [j][i]
            new_face.push_back(edge_point_idx[get_edge(faces_store[j][(i+face_size-1)%face_size], faces_store[j][i])]);
            new_faces.push_back(new_face);
            i++;
        }
        j++;
    }
    clear_mesh();

    // Return the new vertices and faces
    return std::make_pair(new_vertices, new_faces);
}

void mesh::triangulate_mesh(std::vector<std::vector<int>> &faces, int nf){
    // Triangulate the mesh
    for(auto x : faces){
        std::vector<int> temp_face;
        for(auto y : x){
            std::cout << y << ' ';
            temp_face.push_back(y);
        }
        std::cout << '\n';
        faces_store.push_back(temp_face);
    }
    std::vector<std::vector<int>> new_faces;
    int face_idx = 0;
    int tri_face_idx = 0;
    wh(face_idx, nf){
        const std::vector<int>& face = faces[face_idx];
        int face_size = face.size();
        // Boundary Edge dealing
        int j = 0;
        wh(j, face_size){
            boundary_edges.insert(std::make_pair(face[j], face[(j+1)%face_size]));
            // Check if boundary edge is already present in the map
            j++;
        }
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
                face_idx_tri_to_poly[tri_face_idx] = face_idx;
                tri_face_idx++;
                i++;
            }
        } else {
            new_faces.push_back(face);
            triangles_ivec.push_back(glm::ivec3(face[0], face[1], face[2]));
            face_idx_tri_to_poly[tri_face_idx] = face_idx;
            tri_face_idx++;
        }
        face_idx++;
    }
    faces = new_faces;

    for(const auto& edge_ : boundary_edges){
        boundary_edges_vec.emplace_back(edge_.first, edge_.second);
    }

}
        
void mesh::face_set_construction(const std::vector<std::vector<int>> &in_faces, int nf){
    // Get the correct faces using consistent ordering
    // std::cout << " Entered face set construction" << std::endl;
    // // Print input
    // for(auto face : in_faces){
    //     for(auto vertex : face){
    //         std::cout << vertex << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // std::cout << " Number of faces: " << nf << std::endl;


    std::vector<std::vector<int>> new_faces = new_consistent_faces(in_faces, nf);
    // std::vector<std::vector<int>> new_faces = in_faces;

    // std::cout << " Done consistent ordering" << std::endl;

    int face_idx = 0;
    wh(face_idx,nf){

        // std::cout << " Face " << face_idx << ": ";

        // Construction work
        int face_size = new_faces[face_idx].size();
        int i = 0;
        int half_edge_idx = half_edge_vector.size();
        wh(i, face_size){

            // std::cout << "Processing edge " << i << std::endl;

            half_edge new_half_edge;
            // Edge from A to B will point B
            new_half_edge.vertex_idx = new_faces[face_idx][(i+1)%face_size];
            new_half_edge.face_idx = face_idx;
            new_half_edge.next_half_edge_idx = half_edge_idx + (i+1)%face_size;

            std::pair<int, int> twin_edge = get_edge(new_faces[face_idx][(i+1)%face_size], new_faces[face_idx][i]);
            std::pair<int, int> curr_edge = get_edge(new_faces[face_idx][i], new_faces[face_idx][(i+1)%face_size]);

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
            update_vertex(new_half_edge.vertex_idx , half_edge_idx + i);
            i++;
        }
        // Update faces
        face new_face;
        new_face.face_idx = face_idx;
        new_face.half_edge_idx = half_edge_idx;
        // Skip normal
        faces.push_back(new_face);
        face_normal_gen(face_idx);
        face_idx++;
    }

    // Call vertex normal update for all vertices
    int v_iter = 0;
    int v_lim = vertices.size();
    wh(v_iter,v_lim){
        vertex_normal_update(v_iter);
        v_iter++;
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
                  << "Face Area = "<< f.area << ", "
                  << "Face Normal = (" << f.face_normal.x << ", " << f.face_normal.y << ", " << f.face_normal.z << ")" << std::endl;
    }

    // Print edge-to-half-edge mapping
    std::cout << "\nEdge-to-Half-Edge Mapping:" << std::endl;
    for (const auto& entry : edge_to_half_edge) {
        const std::pair<int, int>& edge = entry.first;
        int half_edge_idx = entry.second;
        std::cout << "Edge (" << edge.first << ", " << edge.second << ") -> Half-Edge Index = " << half_edge_idx << std::endl;
    }

    // // Print boundary edges and the mapping to faces
    // std::cout << "\nBoundary Edges and Mapping to Faces:" << std::endl;
    // for (const auto& entry : boundary_edge_to_face) {
    //     const std::pair<int, int>& edge = entry.first;
    //     const std::pair<int, int>& face_pair = entry.second;
    //     std::cout << "Boundary Edge (" << edge.first << ", " << edge.second << ") -> Face Indices = (" << face_pair.first << ", " << face_pair.second << ")" << std::endl;
    // }

    // // Print boundary edges vector
    // std::cout << "\nBoundary Edges:" << std::endl;
    // for (size_t i = 0; i < boundary_edges_vec.size(); ++i) {
    //     const glm::ivec2& edge = boundary_edges_vec[i];
    //     std::cout << "Boundary Edge " << i << ": (" << edge.x << ", " << edge.y << ")" << std::endl;
    // }

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

void mesh::clear_mesh(){
    edge_to_half_edge.clear();
    boundary_edges.clear();
    boundary_edge_to_face.clear();
    boundary_edges_vec.clear();
    faces_store.clear();
    face_idx_tri_to_poly.clear();
    vertices_pos.clear();
    vertices_normal.clear();
    triangles_ivec.clear();
    vertices.clear();
    half_edge_vector.clear();
    faces.clear();
}