#include "viewer.hpp"
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>
#include <chrono>
#include <random>
std::random_device rd;
std::mt19937 gen(rd());

void split_string(std::string s, std::vector<std::string> &res){
    res.clear();
    std::string my = "";
    int i=0;
    while(i < s.size()){
        if(s[i] == ' ' && my.size() > 0){
            res.push_back(my);
            my = "";
        } else if(s[i] != ' '){
            my += s[i];
        }
        i++;
    }
    if(my.size() > 0){
        res.push_back(my);
    }
}

namespace COL781 {
    namespace Viewer {

        namespace GL = COL781::OpenGL;

        void Camera::initialize(float aspect) {
            firstMouse = true;
            yaw   = -90.0f;    
            pitch =  0.0f;
            lastX =  800.0f / 2.0;
            lastY =  600.0 / 2.0;
            fov   =  60.0f;

            this->aspect = aspect;

            position = glm::vec3(0.0f, 0.0f,  1.5f);
            lookAt = glm::vec3(0.0f, 0.0f, 0.0f);
            up = glm::vec3(0.0f, 1.0f,  0.0f);

            updateViewMatrix();
        }

        glm::mat4 Camera::getViewMatrix() {
            return viewMatrix;
        }

        void Camera::updateViewMatrix() {
            viewMatrix = glm::lookAt(position, lookAt, up);
        }

        glm::mat4 Camera::getProjectionMatrix() {
            return glm::perspective(glm::radians(fov), aspect, 0.1f, 100.0f);
        }
            glm::vec3 getRightVector();

        glm::vec3 Camera:: getViewDir() {
            return -glm::transpose(viewMatrix)[2];
        }

        glm::vec3 Camera::getRightVector() {
            return glm::transpose(viewMatrix)[0];
        }

        void Camera::setCameraView(glm::vec3 position_vector, glm::vec3 lookat_vector, glm::vec3 up_vector) {
            position = std::move(position_vector);
            lookAt = std::move(lookat_vector);
            up = std::move(up_vector);

            viewMatrix = glm::lookAt(position, lookAt, up);
        }

        bool Viewer::initialize(const std::string &title, int width, int height) {
            if (!r.initialize(title.c_str(), width, height))
                return false;
            program = r.createShaderProgram(
                r.vsBlinnPhong(),
                r.fsBlinnPhong()
            );
            r.useShaderProgram(program);
            object = r.createObject();
            wireframe = r.createObject();
            r.enableDepthTest();
            camera.initialize((float)width/(float)height);
            return true;
        }

        glm::mat4 calculateModelMatrix(const glm::vec3* vertices, int nv) {
            if (nv <= 0) {
                return glm::mat4(1.0f);
            }
            glm::vec3 min = vertices[0];
            glm::vec3 max = vertices[0];
            for (int i = 1; i < nv; ++i) {
                min = glm::min(min, vertices[i]);
                max = glm::max(max, vertices[i]);
            }
            glm::vec3 center = (min + max) * 0.5f;
            glm::vec3 size = max - min;
            float max_dimension = glm::max(glm::max(size.x, size.y), size.z);
            if (max_dimension <= 0.0f) {
                max_dimension = 1.0f;
            }
            float scale_factor = 1.0f / max_dimension;
            glm::mat4 translation = glm::translate(glm::mat4(1.0f), -center);
            glm::mat4 scaling = glm::scale(glm::mat4(1.0f), glm::vec3(scale_factor));
            return scaling * translation;
        }

        void Viewer::setMesh(int nv, int nt, int ne, const glm::vec3* vertices, const glm::ivec3* triangles, const glm::ivec2* edges, const glm::vec3* normals) {
            if(normals == nullptr) {
                glm::vec3* normalsz = new glm::vec3[nv];
                for(int i = 0; i < nv; i++) {
                    normalsz[i] = glm::vec3(0.0, 0.0, 0.0);
                }
                normals = normalsz;
            }
            r.setVertexAttribs(object, 0, nv, vertices);
            r.setVertexAttribs(object, 1, nv, normals);
            r.setTriangleIndices(object, nt, triangles);
            r.setVertexAttribs(wireframe, 0, nv, vertices);
            r.setVertexAttribs(wireframe, 1, nv, normals);
            r.setEdgeIndices(wireframe, ne, edges); // These are the outline edges of the mesh
            stagetransform = calculateModelMatrix(vertices, nv);

            // Checking working of my a2 functions
            // my_mesh.vertex_set_construction(vertices, normals, nv, true);
            // my_mesh.face_set_construction(my_mesh.tri_converter(triangles,nt), nt);
            // my_mesh.print_ds();

        }

        void Viewer::setMesh_testing(int nv, int nt, int ne, const glm::vec3* vertices, const glm::ivec3* triangles, const glm::ivec2* edges, const glm::vec3* normals) {
            if(normals == nullptr) {
                glm::vec3* normalsz = new glm::vec3[nv];
                for(int i = 0; i < nv; i++) {
                    normalsz[i] = glm::vec3(0.0, 0.0, 0.0);
                }
                normals = normalsz;
            }
            // My implementation
            // First set the vertices and normals
            my_mesh.vertex_set_construction(vertices, normals, nv, true);
            // Then convert the faces to triangles
            std::vector<std::vector<int>> poly_faces = my_mesh.tri_converter(triangles, nt);
            // Then set the faces
            my_mesh.face_set_construction(poly_faces, nt);

            // Now update the Rasterizer
            r.setVertexAttribs(object, 0, my_mesh.get_num_vertices() , my_mesh.get_vertices_pos());
            r.setVertexAttribs(object, 1, my_mesh.get_num_vertices(), my_mesh.get_vertices_normal());
            r.setTriangleIndices(object, my_mesh.get_num_faces(), my_mesh.get_triangles());
            r.setVertexAttribs(wireframe, 0, my_mesh.get_num_vertices(), my_mesh.get_vertices_pos());
            r.setVertexAttribs(wireframe, 1, my_mesh.get_num_vertices(), my_mesh.get_vertices_normal());
            r.setEdgeIndices(wireframe, my_mesh.get_num_boundary_edges(), my_mesh.get_boundary_edges());
            stagetransform = calculateModelMatrix(my_mesh.get_vertices_pos(), my_mesh.get_num_vertices());
        }

        void Viewer::umbrella_update_mesh(float delta, int iters){
            // Normals will be present
            // Call umbrella update all
            my_mesh.umbrella_update_all(delta,iters);
            // Now update the Rasterizer
            r.setVertexAttribs(object, 0, my_mesh.get_num_vertices() , my_mesh.get_vertices_pos());
            r.setVertexAttribs(object, 1, my_mesh.get_num_vertices(), my_mesh.get_vertices_normal());
            r.setTriangleIndices(object, my_mesh.get_num_faces(), my_mesh.get_triangles());
            r.setVertexAttribs(wireframe, 0, my_mesh.get_num_vertices(), my_mesh.get_vertices_pos());
            r.setVertexAttribs(wireframe, 1, my_mesh.get_num_vertices(), my_mesh.get_vertices_normal());
            r.setEdgeIndices(wireframe, my_mesh.get_num_boundary_edges(), my_mesh.get_boundary_edges());
            stagetransform = calculateModelMatrix(my_mesh.get_vertices_pos(), my_mesh.get_num_vertices());
        }

        void Viewer::load_obj_file(const std::string &filepath){
            std::ifstream file(filepath);
            std::string line;
            std::vector<glm::vec3> vertices;
            std::vector<glm::vec3> vertice_normal;
            std::vector<std::vector<int>> faces;
            std::vector<std::string> temp;
            std::cout << "Parsing obj file" << std::endl;
            while (std::getline(file, line)) {
                split_string(line, temp);
                if(temp[0] == "v"){
                    // Assume that the vertices are in the format v x y z
                    vertices.push_back(glm::vec3(std::stof(temp[1]), std::stof(temp[2]), std::stof(temp[3])));
                }
                else if(temp[0] == "vn"){
                    // Assume that the normals are in the format vn x y z
                    vertice_normal.push_back(glm::vec3(std::stof(temp[1]), std::stof(temp[2]), std::stof(temp[3])));
                }
                else if(temp[0] == "f"){
                    // Faces are in the format f v1 v2 v3
                    // or f v1/* v2/* v3/*
                    std::vector<int> new_face;
                    for(int i = 1; i < temp.size(); i++){
                        std::string vertex = temp[i];
                        std::string vertex_index = vertex.substr(0, vertex.find("/"));
                        new_face.push_back(std::stoi(vertex_index) - 1);
                    }
                    faces.push_back(new_face);
                }
                else{
                    std::cout << "Ignoring line: " << line << std::endl;
                }
            }
            // Now set the mesh ensuring 0 normals if not present
            if(vertice_normal.size() == 0){
                // Set normals to 0
                int nv = vertices.size();
                int x = 0;
                while(x < nv){
                    vertice_normal.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
                    x++;
                }
            }
            assert(vertices.size() == vertice_normal.size());
            setMesh_new(vertices.size(), vertices.data(), faces, vertice_normal.data());
        }

        void Viewer::catmull_clark(){
            std::cout << "Catmull clark was called " << std::endl;
            // Call the catmull clark subdivision
            std::pair<std::vector<glm::vec3>, std::vector<std::vector<int>>> new_mesh = my_mesh.catmull_clark_subdivision();

            // Now set the mesh with normals as nullptr
            setMesh_new(new_mesh.first.size(), new_mesh.first.data(), new_mesh.second, nullptr);
            // Compute the normals
            my_mesh.vertex_normal_update_mode_all();
        }

        void Viewer::extrude(int idx, float d){
            // Call the extrude face function
            std::pair<std::vector<glm::vec3>, std::vector<std::vector<int>>> new_mesh = my_mesh.extrude_face(idx, d);
            // Now set the mesh with normals as nullptr
            setMesh_new(new_mesh.first.size(), new_mesh.first.data(), new_mesh.second, nullptr);
            // Compute the normals
            my_mesh.vertex_normal_update_mode_all();
        }

        void Viewer::extrude_region(std::vector<int> face_idx, float d){
            // Call the extrude region function
            std::pair<std::vector<glm::vec3>, std::vector<std::vector<int>>> new_mesh = my_mesh.extrude_region(face_idx, d);
            // Now set the mesh with normals as nullptr
            setMesh_new(new_mesh.first.size(), new_mesh.first.data(), new_mesh.second, nullptr);
            // Compute the normals
            my_mesh.vertex_normal_update_mode_all();
        }


        int Viewer::get_closest_face(glm::vec3 p){
            return my_mesh.get_closest_face(p);
        }

        void Viewer::setMesh_new(int nv, const glm::vec3* vertices, const std::vector<std::vector<int>> &poly_faces, const glm::vec3* normals){
            if(normals == nullptr) {
                glm::vec3* normalsz = new glm::vec3[nv];
                for(int i = 0; i < nv; i++) {
                    normalsz[i] = glm::vec3(0.0, 0.0, 0.0);
                }
                normals = normalsz;
            }
            // My implementation
            // First set the vertices and normals
            std::cout << " Entered setMesh_new" << std::endl;
            my_mesh.vertex_set_construction(vertices, normals, nv, true);
            std::cout << " Done setting vertices" << std::endl;
            // Then convert the faces to triangles
            // Need to deal with const part
            // Create a copy of the faces
            int nf = poly_faces.size();
            std::cout << " Calling triangulate mesh" << std::endl;
            std::vector<std::vector<int>> poly_faces_copy = poly_faces;
            my_mesh.triangulate_mesh(poly_faces_copy, nf);
            std::cout << " Done triangulating mesh" << std::endl;
            // Then set the faces
            nf = poly_faces_copy.size();
            std::cout << " Calling face set construction" << std::endl;
            my_mesh.face_set_construction(poly_faces_copy, nf);
            std::cout << " Done face set construction" << std::endl;

            my_mesh.map_maker();

            // Uncomment for internal variables
            // my_mesh.print_ds();

            // Now update the Rasterizer
            r.setVertexAttribs(object, 0, my_mesh.get_num_vertices() , my_mesh.get_vertices_pos());
            r.setVertexAttribs(object, 1, my_mesh.get_num_vertices(), my_mesh.get_vertices_normal());
            r.setTriangleIndices(object, my_mesh.get_num_faces(), my_mesh.get_triangles());
            r.setVertexAttribs(wireframe, 0, my_mesh.get_num_vertices(), my_mesh.get_vertices_pos());
            r.setVertexAttribs(wireframe, 1, my_mesh.get_num_vertices(), my_mesh.get_vertices_normal());
            r.setEdgeIndices(wireframe, my_mesh.get_num_boundary_edges(), my_mesh.get_boundary_edges());
            stagetransform = calculateModelMatrix(my_mesh.get_vertices_pos(), my_mesh.get_num_vertices());
            std::cout << " Done setting mesh" << std::endl;

            my_mesh.print_ds();
        }

        void Viewer::flip_normals(){
            // Fetch the normals
            glm::vec3* normals = my_mesh.get_vertices_normal();
            // Copy the normals
            int nv = my_mesh.get_num_vertices();
            glm::vec3* new_normals = new glm::vec3[nv];
            int i = 0;
            while(i < nv){
                new_normals[i] = -normals[i];
                i++;
            }
            // Now set the mesh with the new normals
            r.setVertexAttribs(object, 0, my_mesh.get_num_vertices() , my_mesh.get_vertices_pos());
            r.setVertexAttribs(object, 1, my_mesh.get_num_vertices(), new_normals);
            r.setTriangleIndices(object, my_mesh.get_num_faces(), my_mesh.get_triangles());
            r.setVertexAttribs(wireframe, 0, my_mesh.get_num_vertices(), my_mesh.get_vertices_pos());
            r.setVertexAttribs(wireframe, 1, my_mesh.get_num_vertices(), new_normals);
            r.setEdgeIndices(wireframe, my_mesh.get_num_boundary_edges(), my_mesh.get_boundary_edges());
            stagetransform = calculateModelMatrix(my_mesh.get_vertices_pos(), my_mesh.get_num_vertices());
            std::cout << " Done flipping normals" << std::endl;
        }

        void Viewer::create_unit_rectangle(int m, int n){
            // Create a unit rectangle with m rows and n columns
            int idx1 = 0;
            int idx2 = 0;
            int nv = (m+1)*(n+1);
            std::vector<glm::vec3> vertices(nv);
            while(idx1 < m+1){
                idx2 = 0;
                while(idx2 < n+1){
                    vertices[idx1*(n+1) + idx2] = glm::vec3((float)idx2/(float)n, (float)idx1/(float)m, 0.0f);
                    idx2++;
                }
                idx1++;
            }
            // Poly faces
            std::vector<std::vector<int>> poly_faces;
            idx1 = 0;
            while(idx1 < m){
                idx2 = 0;
                while(idx2 < n){
                    // Make a square and add to poly_faces
                    std::vector<int> square;
                    square.push_back(idx1*(n+1) + idx2);
                    square.push_back(idx1*(n+1) + idx2 + 1);
                    square.push_back((idx1+1)*(n+1) + idx2 + 1);
                    square.push_back((idx1+1)*(n+1) + idx2);
                    poly_faces.push_back(square);
                    idx2++;
                }
                idx1++;
            }
            // Now create normals all pointing in 0,0,1 direction
            glm::vec3* normals = new glm::vec3[nv];
            idx1 = 0;
            while(idx1 < nv){
                normals[idx1] = glm::vec3(0.0f, 0.0f, 1.0f);
                idx1++;
            }
            setMesh_new(nv, vertices.data(), poly_faces, normals);
        }

        void Viewer::create_sphere(int slices, int stacks){
            // We will have a total of 
            std::cout << " Entered create sphere" << std::endl;
            int nv = 2 + slices*(stacks-1);
            // Now get the vertices
            std::vector<glm::vec3> vertices;
            std::vector<glm::vec3> normals;
            // Add the top and bottom vertices
            float radius = 0.5f;
            vertices.push_back(glm::vec3(0.0f, 0.0f, radius));
            normals.push_back(glm::vec3(0.0f, 0.0f, 1.0f));
            vertices.push_back(glm::vec3(0.0f, 0.0f, -radius));
            normals.push_back(glm::vec3(0.0f, 0.0f, -1.0f));
            // Now add the vertices in the middle
            float theta = 0.0f;
            float phi = 0.0f;
            float d_theta = 3.14f/(float)stacks;
            float d_phi = 2.0f*3.14f/(float)slices;
            int idx1 = 0;
            int idx2 = 0;
            while(idx1 < stacks-1){
                theta += d_theta;
                phi = 0.0f;
                idx2 = 0;
                while(idx2 < slices){
                    vertices.push_back(glm::vec3(radius*glm::sin(theta)*glm::cos(phi), radius*glm::sin(theta)*glm::sin(phi), radius*glm::cos(theta)));
                    normals.push_back(glm::vec3(glm::sin(theta)*glm::cos(phi), glm::sin(theta)*glm::sin(phi), glm::cos(theta)));
                    phi += d_phi;
                    idx2++;
                }
                idx1++;
            }
            // Now create the faces
            std::vector<std::vector<int>> poly_faces;
            // Pole faces
            // North
            idx1 = 0;
            while(idx1 < slices){
                std::vector<int> pole_face;
                pole_face.push_back(0);
                pole_face.push_back(2 + idx1);
                pole_face.push_back(2 + (idx1+1)%slices);
                poly_faces.push_back(pole_face);
                idx1++;
            }
            // South
            idx1 = 0;
            while(idx1 < slices){
                std::vector<int> pole_face;
                pole_face.push_back(1);
                pole_face.push_back(2 + (stacks-2)*slices + idx1);
                pole_face.push_back(2 + (stacks-2)*slices + (idx1+1)%slices);
                poly_faces.push_back(pole_face);
                idx1++;
            }
            // Now the middle faces
            idx1 = 0;
            while(idx1 < stacks-2){
                idx2 = 0;
                while(idx2 < slices){
                    // Make a square and add to poly_faces
                    std::vector<int> square;
                    square.push_back(2 + idx1*slices + idx2);
                    square.push_back(2 + idx1*slices + (idx2+1)%slices);
                    square.push_back(2 + (idx1+1)*slices + (idx2+1)%slices);
                    square.push_back(2 + (idx1+1)*slices + idx2);
                    poly_faces.push_back(square);
                    idx2++;
                }
                idx1++;
            }

            std::cout << " Called setMesh_new" << std::endl;
            // Set Mesh
            setMesh_new(nv, vertices.data(), poly_faces, normals.data());

        }

        void Viewer::create_cube(int m, int n, int o){
            int nv = 0;

            std::vector<glm::vec3> vertices;
            std::vector<glm::vec3> normals;
            std::vector<std::vector<int>> poly_faces;
            float cube_side = 1.0f;

            // Create (m+1)*(n+1)*(o+1) vertices

            int idx1 = 0;
            int idx2 = 0;
            int idx3 = 0;

            while(idx1 < m+1){
                idx2 = 0;
                while(idx2 < n+1){
                    idx3 = 0;
                    while(idx3 < o+1){
                        vertices.push_back(glm::vec3((float)idx2/(float)n - 0.5f , (float)idx1/(float)m - 0.5f , (float)idx3/(float)o - 0.5f ));
                        normals.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
                        idx3++;
                    }
                    idx2++;
                }
                idx1++;
            }

            // Create faces
            // Front and back faces
            idx1 = 0;
            while(idx1 < m){
                idx2 = 0;
                while(idx2 < n){
                    std::vector<int> front_face;
                    front_face.push_back(idx1*(n+1)*(o+1) + idx2*(o+1));
                    front_face.push_back((idx1+1)*(n+1)*(o+1) + idx2*(o+1));
                    front_face.push_back((idx1+1)*(n+1)*(o+1) + (idx2+1)*(o+1));
                    front_face.push_back(idx1*(n+1)*(o+1) + (idx2+1)*(o+1));

                    std::vector<int> back_face;
                    back_face.push_back(idx1*(n+1)*(o+1) + (idx2+1)*(o+1) + o);
                    back_face.push_back((idx1+1)*(n+1)*(o+1) + (idx2+1)*(o+1) + o);
                    back_face.push_back((idx1+1)*(n+1)*(o+1) + idx2*(o+1) + o);
                    back_face.push_back(idx1*(n+1)*(o+1) + idx2*(o+1) + o);

                    poly_faces.push_back(front_face);
                    poly_faces.push_back(back_face);
                    idx2++;
                }
                idx1++;
            }

            // Top and bottom faces
            idx1 = 0;
            while(idx1 < m){
                idx2 = 0;
                while(idx2 < o){
                    std::vector<int> top_face;
                    top_face.push_back(idx1*(n+1)*(o+1) + n*(o+1) + idx2);
                    top_face.push_back((idx1+1)*(n+1)*(o+1) + n*(o+1) + idx2);
                    top_face.push_back((idx1+1)*(n+1)*(o+1) + n*(o+1) + (idx2+1));
                    top_face.push_back(idx1*(n+1)*(o+1) + n*(o+1) + (idx2+1));

                    std::vector<int> bottom_face;
                    bottom_face.push_back(idx1*(n+1)*(o+1) + (idx2+1));
                    bottom_face.push_back((idx1+1)*(n+1)*(o+1) + (idx2+1));
                    bottom_face.push_back((idx1+1)*(n+1)*(o+1) + idx2);
                    bottom_face.push_back(idx1*(n+1)*(o+1) + idx2);

                    poly_faces.push_back(top_face);
                    poly_faces.push_back(bottom_face);
                    idx2++;
                }
                idx1++;
            }

            // Left and right faces
            idx1 = 0;
            while(idx1 < n){
                idx2 = 0;
                while(idx2 < o){
                    std::vector<int> left_face;
                    left_face.push_back(idx1*(o+1) + idx2);
                    left_face.push_back((idx1+1)*(o+1) + idx2);
                    left_face.push_back((idx1+1)*(o+1) + (idx2+1));
                    left_face.push_back(idx1*(o+1) + (idx2+1));

                    std::vector<int> right_face;
                    right_face.push_back(m*(n+1)*(o+1) + idx1*(o+1) + idx2);
                    right_face.push_back(m*(n+1)*(o+1) + (idx1+1)*(o+1) + idx2);
                    right_face.push_back(m*(n+1)*(o+1) + (idx1+1)*(o+1) + (idx2+1));
                    right_face.push_back(m*(n+1)*(o+1) + idx1*(o+1) + (idx2+1));

                    poly_faces.push_back(left_face);
                    poly_faces.push_back(right_face);
                    idx2++;
                }
                idx1++;
            }

            // Correcting normals
            // Front and back faces have normals in z direction
            idx1 = 0;
            while(idx1 < m+1){
                idx2 = 0;
                while(idx2 < n+1){
                    // Correction
                    normals[idx1*(n+1)*(o+1) + idx2*(o+1)] = glm::vec3(0.0f, 0.0f, -1.0f);
                    normals[idx1*(n+1)*(o+1) + idx2*(o+1) + o] = glm::vec3(0.0f, 0.0f, 1.0f);
                    idx2++;
                }
                idx1++;
            }

            // Top and bottom faces have normals in y direction
            idx1 = 0;
            while(idx1 < m+1){
                idx2 = 0;
                while(idx2 < o+1){
                    normals[idx1*(n+1)*(o+1) + n*(o+1) + idx2] = glm::vec3(0.0f, 1.0f, 0.0f);
                    normals[idx1*(n+1)*(o+1) + idx2] = glm::vec3(0.0f, -1.0f, 0.0f);
                    idx2++;
                }
                idx1++;
            }

            // Left and right faces have normals in x direction
            idx1 = 0;
            while(idx1 < n+1){
                idx2 = 0;
                while(idx2 < o+1){
                    normals[idx1*(o+1) + idx2] = glm::vec3(-1.0f, 0.0f, 0.0f);
                    normals[m*(n+1)*(o+1) + idx1*(o+1) + idx2] = glm::vec3(1.0f, 0.0f, 0.0f);
                    idx2++;
                }
                idx1++;
            }
            // Edge vertice normals are between the two faces
            // x - direction edges
            idx1 = 0;
            while(idx1 < m+1){
                // for the 4 edges
                normals[idx1*(n+1)*(o+1)] = glm::normalize(glm::vec3(0.0f, -1.0f, 1.0f));
                normals[idx1*(n+1)*(o+1) + o] = glm::normalize(glm::vec3(0.0f, -1.0f, -1.0f));
                normals[idx1*(n+1)*(o+1) + n*(o+1)] = glm::normalize(glm::vec3(0.0f, 1.0f, 1.0f));
                normals[idx1*(n+1)*(o+1) + n*(o+1) + o] = glm::normalize(glm::vec3(0.0f, 1.0f, -1.0f));
                idx1++; 
            }
            // y - direction edges
            idx1 = 0;
            while(idx1 < n+1){
                // for the 4 edges
                normals[idx1*(o+1)] = glm::normalize(glm::vec3(-1.0f, 0.0f, -1.0f));
                normals[idx1*(o+1) + m*(n+1)*(o+1)] = glm::normalize(glm::vec3(1.0f, 0.0f, -1.0f));
                normals[idx1*(o+1) + o] = glm::normalize(glm::vec3(-1.0f, 0.0f, 1.0f));
                normals[idx1*(o+1) + m*(n+1)*(o+1) + o] = glm::normalize(glm::vec3(1.0f, 0.0f, 1.0f));
                idx1++; 
            }
            // z - direction edges
            idx1 = 0;
            while(idx1 < o+1){
                // for the 4 edges
                normals[idx1] = glm::normalize(glm::vec3(-1.0f, -1.0f, 0.0f));
                normals[idx1 + m*(n+1)*(o+1)] = glm::normalize(glm::vec3(1.0f, -1.0f, 0.0f));
                normals[idx1 + n*(o+1)] = glm::normalize(glm::vec3(-1.0f, 1.0f, 0.0f));
                normals[idx1 + m*(n+1)*(o+1) + n*(o+1)] = glm::normalize(glm::vec3(1.0f, 1.0f, 0.0f));
                idx1++; 
            }

            // Corner normals
            normals[0] = glm::normalize(glm::vec3(-1.0f, -1.0f, -1.0f));
            normals[m*(n+1)*(o+1)] = glm::normalize(glm::vec3(1.0f, -1.0f, -1.0f));
            normals[n*(o+1)] = glm::normalize(glm::vec3(-1.0f, 1.0f, -1.0f));
            normals[m*(n+1)*(o+1) + n*(o+1)] = glm::normalize(glm::vec3(1.0f, 1.0f, -1.0f));
            normals[o] = glm::normalize(glm::vec3(-1.0f, -1.0f, 1.0f));
            normals[m*(n+1)*(o+1) + o] = glm::normalize(glm::vec3(1.0f, -1.0f, 1.0f));
            normals[n*(o+1) + o] = glm::normalize(glm::vec3(-1.0f, 1.0f, 1.0f));
            normals[m*(n+1)*(o+1) + n*(o+1) + o] = glm::normalize(glm::vec3(1.0f, 1.0f, 1.0f));

            // Create a set of used vertices
            std::set<int> used_vertices;
            for(const auto& face : poly_faces){
                for(const auto& vertex : face){
                    used_vertices.insert(vertex);
                }
            }

            nv = 0;
            // Create a map from old vertices to new vertices
            std::map<int, int> old_to_new;
            for(const auto& vertex : used_vertices){
                old_to_new[vertex] = nv;
                nv++;
            }

            // Now create the new vertices and normals
            std::vector<glm::vec3> new_vertices(nv);
            std::vector<glm::vec3> new_normals(nv);
            for(const auto& vertex : used_vertices){
                new_vertices[old_to_new[vertex]] = vertices[vertex];
                new_normals[old_to_new[vertex]] = normals[vertex];
            }

            // Now update the faces
            for(auto& face : poly_faces){
                for(auto& vertex : face){
                    vertex = old_to_new[vertex];
                }
            }

            setMesh_new(nv, new_vertices.data(), poly_faces, new_normals.data());
            my_mesh.vertex_normal_update_mode_all();

        }
        
        void Viewer::create_cube_new(int m, int n, int o){
            // Create (m+1)*(n+1)*(o+1) vertices
            // Centered at origin
            int idx1 = 0; // for x
            int idx2 = 0; // for y
            int idx3 = 0; // for z
        
            // Create vertices
            std::vector<glm::vec3> vertices;
            while(idx1 < m+1){
                idx2 = 0;
                while(idx2 < n+1){
                    idx3 = 0;
                    while(idx3 < o+1){
                        vertices.push_back(glm::vec3(
                            (float)idx1/(float)m - 0.5f, 
                            (float)idx2/(float)n - 0.5f, 
                            (float)idx3/(float)o - 0.5f
                        ));
                        idx3++;
                    }
                    idx2++;
                }
                idx1++;
            }
        
            // Create faces
            std::vector<std::vector<int>> poly_faces;
        
            // Front face (Positive z)
            idx2 = 0;
            while(idx2 < n){
                idx1 = 0;
                while(idx1 < m){
                    std::vector<int> front_face;
                    int base = idx1 + idx2 * (m+1) + o * (m+1) * (n+1);  // z = o
                    front_face.push_back(base);                    // bottom-left
                    front_face.push_back(base + (m+1));           // bottom-right
                    front_face.push_back(base + (m+1) + 1);       // top-right
                    front_face.push_back(base + 1);               // top-left
                    poly_faces.push_back(front_face);
                    idx1++;
                }
                idx2++;
            }
        
            // Back face (Negative z)
            idx2 = 0;
            while(idx2 < n){
                idx1 = 0;
                while(idx1 < m){
                    std::vector<int> back_face;
                    int base = idx1 + idx2 * (m+1);  // z = 0
                    back_face.push_back(base);                    // bottom-left
                    back_face.push_back(base + 1);                // bottom-right
                    back_face.push_back(base + (m+1) + 1);        // top-right
                    back_face.push_back(base + (m+1));            // top-left
                    poly_faces.push_back(back_face);
                    idx1++;
                }
                idx2++;
            }
        
            // Left face (Negative x)
            idx2 = 0;
            while(idx2 < n){
                idx3 = 0;
                while(idx3 < o){
                    std::vector<int> left_face;
                    int base = idx2 * (m+1) + idx3 * (m+1) * (n+1);  // x = 0
                    left_face.push_back(base);                    // bottom-left
                    left_face.push_back(base + (m+1)*(n+1));      // bottom-right
                    left_face.push_back(base + (m+1)*(n+1) + (m+1));  // top-right
                    left_face.push_back(base + (m+1));            // top-left
                    poly_faces.push_back(left_face);
                    idx3++;
                }
                idx2++;
            }
        
            // Right face (Positive x)
            idx2 = 0;
            while(idx2 < n){
                idx3 = 0;
                while(idx3 < o){
                    std::vector<int> right_face;
                    int base = m + idx2 * (m+1) + idx3 * (m+1) * (n+1);  // x = m
                    right_face.push_back(base);                   // bottom-left
                    right_face.push_back(base + (m+1));           // bottom-right
                    right_face.push_back(base + (m+1)*(n+1) + (m+1));  // top-right
                    right_face.push_back(base + (m+1)*(n+1));     // top-left
                    poly_faces.push_back(right_face);
                    idx3++;
                }
                idx2++;
            }
        
            // Bottom face (Negative y)
            idx1 = 0;
            while(idx1 < m){
                idx3 = 0;
                while(idx3 < o){
                    std::vector<int> bottom_face;
                    int base = idx1 + idx3 * (m+1) * (n+1);  // y = 0
                    bottom_face.push_back(base);                  // bottom-left
                    bottom_face.push_back(base + (m+1)*(n+1));    // bottom-right
                    bottom_face.push_back(base + (m+1)*(n+1) + 1); // top-right
                    bottom_face.push_back(base + 1);              // top-left
                    poly_faces.push_back(bottom_face);
                    idx3++;
                }
                idx1++;
            }
        
            // Top face (Positive y)
            idx1 = 0;
            while(idx1 < m){
                idx3 = 0;
                while(idx3 < o){
                    std::vector<int> top_face;
                    int base = idx1 + n * (m+1) + idx3 * (m+1) * (n+1);  // y = n
                    top_face.push_back(base);                    // bottom-left
                    top_face.push_back(base + 1);                // bottom-right
                    top_face.push_back(base + (m+1)*(n+1) + 1);   // top-right
                    top_face.push_back(base + (m+1)*(n+1));      // top-left
                    poly_faces.push_back(top_face);
                    idx3++;
                }
                idx1++;
            }

            // Create a set of used vertices
            std::set<int> used_vertices;
            for(const auto& face : poly_faces){
                for(const auto& vertex : face){
                    used_vertices.insert(vertex);
                }
            }

            int nv = 0;
            // Create a map from old vertices to new vertices
            std::map<int, int> old_to_new;
            for(const auto& vertex : used_vertices){
                old_to_new[vertex] = nv;
                nv++;
            }

            // Now create the new vertices and new faces
            std::vector<glm::vec3> new_vertices(nv);
            std::vector<glm::vec3> new_normals(nv);
            for(const auto& vertex : used_vertices){
                new_vertices[old_to_new[vertex]] = vertices[vertex];
                // Setting the normals 
                if(vertices[vertex].x == -0.5f){
                    new_normals[old_to_new[vertex]] = glm::vec3(-1.0f, 0.0f, 0.0f);
                }
                else if(vertices[vertex].x == 0.5f){
                    new_normals[old_to_new[vertex]] = glm::vec3(1.0f, 0.0f, 0.0f);
                }
                else if(vertices[vertex].y == -0.5f){
                    new_normals[old_to_new[vertex]] = glm::vec3(0.0f, -1.0f, 0.0f);
                }
                else if(vertices[vertex].y == 0.5f){
                    new_normals[old_to_new[vertex]] = glm::vec3(0.0f, 1.0f, 0.0f);
                }
                else if(vertices[vertex].z == -0.5f){
                    new_normals[old_to_new[vertex]] = glm::vec3(0.0f, 0.0f, -1.0f);
                }
                else if(vertices[vertex].z == 0.5f){
                    new_normals[old_to_new[vertex]] = glm::vec3(0.0f, 0.0f, 1.0f);
                }
                else{
                    new_normals[old_to_new[vertex]] = glm::vec3(0.0f, 0.0f, 0.0f);
                }

                // For the edges take the average of the normals of the two faces
                if(vertices[vertex].x == 0.5f && vertices[vertex].y == 0.5f){
                    new_normals[old_to_new[vertex]] = glm::normalize(new_normals[old_to_new[vertex]] + glm::vec3(0.0f, 1.0f, 1.0f));
                }
                else if(vertices[vertex].x == -0.5f && vertices[vertex].y == 0.5f){
                    new_normals[old_to_new[vertex]] = glm::normalize(new_normals[old_to_new[vertex]] + glm::vec3(0.0f, 1.0f, -1.0f));
                }
                else if(vertices[vertex].x == 0.5f && vertices[vertex].y == -0.5f){
                    new_normals[old_to_new[vertex]] = glm::normalize(new_normals[old_to_new[vertex]] + glm::vec3(0.0f, -1.0f, 1.0f));
                }
                else if(vertices[vertex].x == -0.5f && vertices[vertex].y == -0.5f){
                    new_normals[old_to_new[vertex]] = glm::normalize(new_normals[old_to_new[vertex]] + glm::vec3(0.0f, -1.0f, -1.0f));
                }
                else if(vertices[vertex].x == 0.5f && vertices[vertex].z == 0.5f){
                    new_normals[old_to_new[vertex]] = glm::normalize(new_normals[old_to_new[vertex]] + glm::vec3(1.0f, 0.0f, 1.0f));
                }
                else if(vertices[vertex].x == -0.5f && vertices[vertex].z == 0.5f){
                    new_normals[old_to_new[vertex]] = glm::normalize(new_normals[old_to_new[vertex]] + glm::vec3(-1.0f, 0.0f, 1.0f));
                }
                else if(vertices[vertex].x == 0.5f && vertices[vertex].z == -0.5f){
                    new_normals[old_to_new[vertex]] = glm::normalize(new_normals[old_to_new[vertex]] + glm::vec3(1.0f, 0.0f, -1.0f));
                }
                else if(vertices[vertex].x == -0.5f && vertices[vertex].z == -0.5f){
                    new_normals[old_to_new[vertex]] = glm::normalize(new_normals[old_to_new[vertex]] + glm::vec3(-1.0f, 0.0f, -1.0f));
                }
                else if(vertices[vertex].y == 0.5f && vertices[vertex].z == 0.5f){
                    new_normals[old_to_new[vertex]] = glm::normalize(new_normals[old_to_new[vertex]] + glm::vec3(0.0f, 1.0f, 1.0f));
                }
                else if(vertices[vertex].y == -0.5f && vertices[vertex].z == 0.5f){
                    new_normals[old_to_new[vertex]] = glm::normalize(new_normals[old_to_new[vertex]] + glm::vec3(0.0f, -1.0f, 1.0f));
                }
                else if(vertices[vertex].y == 0.5f && vertices[vertex].z == -0.5f){
                    new_normals[old_to_new[vertex]] = glm::normalize(new_normals[old_to_new[vertex]] + glm::vec3(0.0f, 1.0f, -1.0f));
                }
                else if(vertices[vertex].y == -0.5f && vertices[vertex].z == -0.5f){
                    new_normals[old_to_new[vertex]] = glm::normalize(new_normals[old_to_new[vertex]] + glm::vec3(0.0f, -1.0f, -1.0f));
                }

                // For the corners take the average of the normals of the three faces
                if(vertices[vertex].x == 0.5f && vertices[vertex].y == 0.5f && vertices[vertex].z == 0.5f){
                    new_normals[old_to_new[vertex]] = glm::normalize(new_normals[old_to_new[vertex]] + glm::vec3(1.0f, 1.0f, 1.0f));
                }
                else if(vertices[vertex].x == -0.5f && vertices[vertex].y == 0.5f && vertices[vertex].z == 0.5f){
                    new_normals[old_to_new[vertex]] = glm::normalize(new_normals[old_to_new[vertex]] + glm::vec3(-1.0f, 1.0f, 1.0f));
                }
                else if(vertices[vertex].x == 0.5f && vertices[vertex].y == -0.5f && vertices[vertex].z == 0.5f){
                    new_normals[old_to_new[vertex]] = glm::normalize(new_normals[old_to_new[vertex]] + glm::vec3(1.0f, -1.0f, 1.0f));
                }
                else if(vertices[vertex].x == -0.5f && vertices[vertex].y == -0.5f && vertices[vertex].z == 0.5f){
                    new_normals[old_to_new[vertex]] = glm::normalize(new_normals[old_to_new[vertex]] + glm::vec3(-1.0f, -1.0f, 1.0f));
                }
                else if(vertices[vertex].x == 0.5f && vertices[vertex].y == 0.5f && vertices[vertex].z == -0.5f){
                    new_normals[old_to_new[vertex]] = glm::normalize(new_normals[old_to_new[vertex]] + glm::vec3(1.0f, 1.0f, -1.0f));
                }
                else if(vertices[vertex].x == -0.5f && vertices[vertex].y == 0.5f && vertices[vertex].z == -0.5f){
                    new_normals[old_to_new[vertex]] = glm::normalize(new_normals[old_to_new[vertex]] + glm::vec3(-1.0f, 1.0f, -1.0f));
                }
                else if(vertices[vertex].x == 0.5f && vertices[vertex].y == -0.5f && vertices[vertex].z == -0.5f){
                    new_normals[old_to_new[vertex]] = glm::normalize(new_normals[old_to_new[vertex]] + glm::vec3(1.0f, -1.0f, -1.0f));
                }
                else if(vertices[vertex].x == -0.5f && vertices[vertex].y == -0.5f && vertices[vertex].z == -0.5f){
                    new_normals[old_to_new[vertex]] = glm::normalize(new_normals[old_to_new[vertex]] + glm::vec3(-1.0f, -1.0f, -1.0f));
                }


            }
            // Now update the faces
            for(auto& face : poly_faces){
                for(auto& vertex : face){
                    vertex = old_to_new[vertex];
                }
            }

            setMesh_new(nv, new_vertices.data(), poly_faces, new_normals.data());
        }
        
        void Viewer::create_noisy_cube(int m, int n, int o){
            int nv = 0;

            std::vector<glm::vec3> vertices;
            std::vector<glm::vec3> normals;
            std::vector<std::vector<int>> poly_faces;
            float cube_side = 1.0f;

            // Create (m+1)*(n+1)*(o+1) vertices
            float noise_lim  = 0.4f;
            std::uniform_real_distribution<float> dist_m(-noise_lim/(float)m, noise_lim/(float)m);
            std::uniform_real_distribution<float> dist_n(-noise_lim/(float)n, noise_lim/(float)n);
            std::uniform_real_distribution<float> dist_o(-noise_lim/(float)o, noise_lim/(float)o);

            int idx1 = 0;
            int idx2 = 0;
            int idx3 = 0;

            glm::vec3 noise;

            while(idx1 < m+1){
                idx2 = 0;
                while(idx2 < n+1){
                    idx3 = 0;
                    while(idx3 < o+1){
                        noise = glm::vec3(dist_m(gen),dist_n(gen),dist_o(gen));
                        vertices.push_back(glm::vec3((float)idx2/(float)n, (float)idx1/(float)m, (float)idx3/(float)o) + noise);
                        normals.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
                        idx3++;
                    }
                    idx2++;
                }
                idx1++;
            }

            // Create faces
            // Front and back faces
            idx1 = 0;
            while(idx1 < m){
                idx2 = 0;
                while(idx2 < n){
                    std::vector<int> front_face;
                    front_face.push_back(idx1*(n+1)*(o+1) + idx2*(o+1));
                    front_face.push_back((idx1+1)*(n+1)*(o+1) + idx2*(o+1));
                    front_face.push_back((idx1+1)*(n+1)*(o+1) + (idx2+1)*(o+1));
                    front_face.push_back(idx1*(n+1)*(o+1) + (idx2+1)*(o+1));

                    std::vector<int> back_face;
                    back_face.push_back(idx1*(n+1)*(o+1) + (idx2+1)*(o+1) + o);
                    back_face.push_back((idx1+1)*(n+1)*(o+1) + (idx2+1)*(o+1) + o);
                    back_face.push_back((idx1+1)*(n+1)*(o+1) + idx2*(o+1) + o);
                    back_face.push_back(idx1*(n+1)*(o+1) + idx2*(o+1) + o);

                    poly_faces.push_back(front_face);
                    poly_faces.push_back(back_face);
                    idx2++;
                }
                idx1++;
            }

            // Top and bottom faces
            idx1 = 0;
            while(idx1 < m){
                idx2 = 0;
                while(idx2 < o){
                    std::vector<int> top_face;
                    top_face.push_back(idx1*(n+1)*(o+1) + n*(o+1) + idx2);
                    top_face.push_back((idx1+1)*(n+1)*(o+1) + n*(o+1) + idx2);
                    top_face.push_back((idx1+1)*(n+1)*(o+1) + n*(o+1) + (idx2+1));
                    top_face.push_back(idx1*(n+1)*(o+1) + n*(o+1) + (idx2+1));

                    std::vector<int> bottom_face;
                    bottom_face.push_back(idx1*(n+1)*(o+1) + (idx2+1));
                    bottom_face.push_back((idx1+1)*(n+1)*(o+1) + (idx2+1));
                    bottom_face.push_back((idx1+1)*(n+1)*(o+1) + idx2);
                    bottom_face.push_back(idx1*(n+1)*(o+1) + idx2);

                    poly_faces.push_back(top_face);
                    poly_faces.push_back(bottom_face);
                    idx2++;
                }
                idx1++;
            }

            // Left and right faces
            idx1 = 0;
            while(idx1 < n){
                idx2 = 0;
                while(idx2 < o){
                    std::vector<int> left_face;
                    left_face.push_back(idx1*(o+1) + idx2);
                    left_face.push_back((idx1+1)*(o+1) + idx2);
                    left_face.push_back((idx1+1)*(o+1) + (idx2+1));
                    left_face.push_back(idx1*(o+1) + (idx2+1));

                    std::vector<int> right_face;
                    right_face.push_back(m*(n+1)*(o+1) + idx1*(o+1) + idx2);
                    right_face.push_back(m*(n+1)*(o+1) + (idx1+1)*(o+1) + idx2);
                    right_face.push_back(m*(n+1)*(o+1) + (idx1+1)*(o+1) + (idx2+1));
                    right_face.push_back(m*(n+1)*(o+1) + idx1*(o+1) + (idx2+1));

                    poly_faces.push_back(left_face);
                    poly_faces.push_back(right_face);
                    idx2++;
                }
                idx1++;
            }

            // Correcting normals
            // Front and back faces have normals in z direction
            idx1 = 0;
            while(idx1 < m+1){
                idx2 = 0;
                while(idx2 < n+1){
                    normals[idx1*(n+1)*(o+1) + idx2*(o+1)] = glm::vec3(0.0f, 0.0f, 1.0f);
                    normals[idx1*(n+1)*(o+1) + idx2*(o+1) + o] = glm::vec3(0.0f, 0.0f, -1.0f);
                    idx2++;
                }
                idx1++;
            }

            // Top and bottom faces have normals in y direction
            idx1 = 0;
            while(idx1 < m+1){
                idx2 = 0;
                while(idx2 < o+1){
                    normals[idx1*(n+1)*(o+1) + n*(o+1) + idx2] = glm::vec3(0.0f, 1.0f, 0.0f);
                    normals[idx1*(n+1)*(o+1) + idx2] = glm::vec3(0.0f, -1.0f, 0.0f);
                    idx2++;
                }
                idx1++;
            }

            // Left and right faces have normals in x direction
            idx1 = 0;
            while(idx1 < n+1){
                idx2 = 0;
                while(idx2 < o+1){
                    normals[idx1*(o+1) + idx2] = glm::vec3(-1.0f, 0.0f, 0.0f);
                    normals[m*(n+1)*(o+1) + idx1*(o+1) + idx2] = glm::vec3(1.0f, 0.0f, 0.0f);
                    idx2++;
                }
                idx1++;
            }
            // Edge vertice normals are between the two faces
            // x - direction edges
            idx1 = 0;
            while(idx1 < m+1){
                // for the 4 edges
                normals[idx1*(n+1)*(o+1)] = glm::normalize(glm::vec3(0.0f, -1.0f, 1.0f));
                normals[idx1*(n+1)*(o+1) + o] = glm::normalize(glm::vec3(0.0f, -1.0f, -1.0f));
                normals[idx1*(n+1)*(o+1) + n*(o+1)] = glm::normalize(glm::vec3(0.0f, 1.0f, 1.0f));
                normals[idx1*(n+1)*(o+1) + n*(o+1) + o] = glm::normalize(glm::vec3(0.0f, 1.0f, -1.0f));
                idx1++; 
            }
            // y - direction edges
            idx1 = 0;
            while(idx1 < n+1){
                // for the 4 edges
                normals[idx1*(o+1)] = glm::normalize(glm::vec3(-1.0f, 0.0f, -1.0f));
                normals[idx1*(o+1) + m*(n+1)*(o+1)] = glm::normalize(glm::vec3(1.0f, 0.0f, -1.0f));
                normals[idx1*(o+1) + o] = glm::normalize(glm::vec3(-1.0f, 0.0f, 1.0f));
                normals[idx1*(o+1) + m*(n+1)*(o+1) + o] = glm::normalize(glm::vec3(1.0f, 0.0f, 1.0f));
                idx1++; 
            }
            // z - direction edges
            idx1 = 0;
            while(idx1 < o+1){
                // for the 4 edges
                normals[idx1] = glm::normalize(glm::vec3(-1.0f, -1.0f, 0.0f));
                normals[idx1 + m*(n+1)*(o+1)] = glm::normalize(glm::vec3(1.0f, -1.0f, 0.0f));
                normals[idx1 + n*(o+1)] = glm::normalize(glm::vec3(-1.0f, 1.0f, 0.0f));
                normals[idx1 + m*(n+1)*(o+1) + n*(o+1)] = glm::normalize(glm::vec3(1.0f, 1.0f, 0.0f));
                idx1++; 
            }

            // Corner normals
            normals[0] = glm::normalize(glm::vec3(-1.0f, -1.0f, -1.0f));
            normals[m*(n+1)*(o+1)] = glm::normalize(glm::vec3(1.0f, -1.0f, -1.0f));
            normals[n*(o+1)] = glm::normalize(glm::vec3(-1.0f, 1.0f, -1.0f));
            normals[m*(n+1)*(o+1) + n*(o+1)] = glm::normalize(glm::vec3(1.0f, 1.0f, -1.0f));
            normals[o] = glm::normalize(glm::vec3(-1.0f, -1.0f, 1.0f));
            normals[m*(n+1)*(o+1) + o] = glm::normalize(glm::vec3(1.0f, -1.0f, 1.0f));
            normals[n*(o+1) + o] = glm::normalize(glm::vec3(-1.0f, 1.0f, 1.0f));
            normals[m*(n+1)*(o+1) + n*(o+1) + o] = glm::normalize(glm::vec3(1.0f, 1.0f, 1.0f));

            // Create a set of used vertices
            std::set<int> used_vertices;
            for(const auto& face : poly_faces){
                for(const auto& vertex : face){
                    used_vertices.insert(vertex);
                }
            }

            nv = 0;
            // Create a map from old vertices to new vertices
            std::map<int, int> old_to_new;
            for(const auto& vertex : used_vertices){
                old_to_new[vertex] = nv;
                nv++;
            }

            // Now create the new vertices and normals
            std::vector<glm::vec3> new_vertices(nv);
            std::vector<glm::vec3> new_normals(nv);
            for(const auto& vertex : used_vertices){
                new_vertices[old_to_new[vertex]] = vertices[vertex];
                new_normals[old_to_new[vertex]] = normals[vertex];
            }

            // Now update the faces
            for(auto& face : poly_faces){
                for(auto& vertex : face){
                    vertex = old_to_new[vertex];
                }
            }

            setMesh_new(nv, new_vertices.data(), poly_faces, new_normals.data());

        }
        
        void Viewer::create_noisy_cube_new(int m, int n, int o){
            // Create (m+1)*(n+1)*(o+1) vertices
            // Centered at origin
            int idx1 = 0; // for x
            int idx2 = 0; // for y
            int idx3 = 0; // for z
            float cube_side = 1.0f;

            // Create (m+1)*(n+1)*(o+1) vertices
            float noise_lim  = 0.2f;
            std::uniform_real_distribution<float> dist_m(-noise_lim/(float)m, noise_lim/(float)m);
            std::uniform_real_distribution<float> dist_n(-noise_lim/(float)n, noise_lim/(float)n);
            std::uniform_real_distribution<float> dist_o(-noise_lim/(float)o, noise_lim/(float)o);
            // Create vertices
            std::vector<glm::vec3> vertices;
            glm::vec3 noise;
            while(idx1 < m+1){
                idx2 = 0;
                while(idx2 < n+1){
                    idx3 = 0;
                    while(idx3 < o+1){
                        noise = glm::vec3(dist_m(gen),dist_n(gen),dist_o(gen));
                        vertices.push_back(glm::vec3(
                            (float)idx1/(float)m - 0.5f, 
                            (float)idx2/(float)n - 0.5f, 
                            (float)idx3/(float)o - 0.5f
                        ) + noise);
                        idx3++;
                    }
                    idx2++;
                }
                idx1++;
            }
        
            // Create faces
            std::vector<std::vector<int>> poly_faces;
        
            // Front face (Positive z)
            idx2 = 0;
            while(idx2 < n){
                idx1 = 0;
                while(idx1 < m){
                    std::vector<int> front_face;
                    int base = idx1 + idx2 * (m+1) + o * (m+1) * (n+1);  // z = o
                    front_face.push_back(base);                    // bottom-left
                    front_face.push_back(base + (m+1));           // bottom-right
                    front_face.push_back(base + (m+1) + 1);       // top-right
                    front_face.push_back(base + 1);               // top-left
                    poly_faces.push_back(front_face);
                    idx1++;
                }
                idx2++;
            }
        
            // Back face (Negative z)
            idx2 = 0;
            while(idx2 < n){
                idx1 = 0;
                while(idx1 < m){
                    std::vector<int> back_face;
                    int base = idx1 + idx2 * (m+1);  // z = 0
                    back_face.push_back(base);                    // bottom-left
                    back_face.push_back(base + 1);                // bottom-right
                    back_face.push_back(base + (m+1) + 1);        // top-right
                    back_face.push_back(base + (m+1));            // top-left
                    poly_faces.push_back(back_face);
                    idx1++;
                }
                idx2++;
            }
        
            // Left face (Negative x)
            idx2 = 0;
            while(idx2 < n){
                idx3 = 0;
                while(idx3 < o){
                    std::vector<int> left_face;
                    int base = idx2 * (m+1) + idx3 * (m+1) * (n+1);  // x = 0
                    left_face.push_back(base);                    // bottom-left
                    left_face.push_back(base + (m+1)*(n+1));      // bottom-right
                    left_face.push_back(base + (m+1)*(n+1) + (m+1));  // top-right
                    left_face.push_back(base + (m+1));            // top-left
                    poly_faces.push_back(left_face);
                    idx3++;
                }
                idx2++;
            }
        
            // Right face (Positive x)
            idx2 = 0;
            while(idx2 < n){
                idx3 = 0;
                while(idx3 < o){
                    std::vector<int> right_face;
                    int base = m + idx2 * (m+1) + idx3 * (m+1) * (n+1);  // x = m
                    right_face.push_back(base);                   // bottom-left
                    right_face.push_back(base + (m+1));           // bottom-right
                    right_face.push_back(base + (m+1)*(n+1) + (m+1));  // top-right
                    right_face.push_back(base + (m+1)*(n+1));     // top-left
                    poly_faces.push_back(right_face);
                    idx3++;
                }
                idx2++;
            }
        
            // Bottom face (Negative y)
            idx1 = 0;
            while(idx1 < m){
                idx3 = 0;
                while(idx3 < o){
                    std::vector<int> bottom_face;
                    int base = idx1 + idx3 * (m+1) * (n+1);  // y = 0
                    bottom_face.push_back(base);                  // bottom-left
                    bottom_face.push_back(base + (m+1)*(n+1));    // bottom-right
                    bottom_face.push_back(base + (m+1)*(n+1) + 1); // top-right
                    bottom_face.push_back(base + 1);              // top-left
                    poly_faces.push_back(bottom_face);
                    idx3++;
                }
                idx1++;
            }
        
            // Top face (Positive y)
            idx1 = 0;
            while(idx1 < m){
                idx3 = 0;
                while(idx3 < o){
                    std::vector<int> top_face;
                    int base = idx1 + n * (m+1) + idx3 * (m+1) * (n+1);  // y = n
                    top_face.push_back(base);                    // bottom-left
                    top_face.push_back(base + 1);                // bottom-right
                    top_face.push_back(base + (m+1)*(n+1) + 1);   // top-right
                    top_face.push_back(base + (m+1)*(n+1));      // top-left
                    poly_faces.push_back(top_face);
                    idx3++;
                }
                idx1++;
            }

            // Create a set of used vertices
            std::set<int> used_vertices;
            for(const auto& face : poly_faces){
                for(const auto& vertex : face){
                    used_vertices.insert(vertex);
                }
            }

            int nv = 0;
            // Create a map from old vertices to new vertices
            std::map<int, int> old_to_new;
            for(const auto& vertex : used_vertices){
                old_to_new[vertex] = nv;
                nv++;
            }

            // Now create the new vertices and new faces
            std::vector<glm::vec3> new_vertices(nv);
            std::vector<glm::vec3> new_normals(nv);
            for(const auto& vertex : used_vertices){
                new_vertices[old_to_new[vertex]] = vertices[vertex];
                new_normals[old_to_new[vertex]] = glm::vec3(0.0f, 0.0f, 0.0f);
            }
            // Now update the faces
            for(auto& face : poly_faces){
                for(auto& vertex : face){
                    vertex = old_to_new[vertex];
                }
            }

            setMesh_new(nv, new_vertices.data(), poly_faces, new_normals.data());
        }

        void Viewer::view(bool flag) {
            // The transformation matrix.
            glm::mat4 model = stagetransform;
            // glm::mat4 model = glm::mat4(1.0f);
            glm::mat4 view;    
            glm::mat4 projection = camera.getProjectionMatrix();

            float deltaAngleX = 2.0 * 3.14 / 800.0;
            float deltaAngleY = 3.14 / 600.0;

            int lastxPos, lastyPos, xPos, yPos;

            SDL_GetMouseState(&lastxPos, &lastyPos);

            auto time_now = std::chrono::high_resolution_clock::now();
            const float delta = 0.04f;
            constexpr std::chrono::milliseconds update_interval(1000); 

            while (!r.shouldQuit()) {

                if(flag == true){
                    // Check if we need to update the mesh
                    auto time_elapsed = std::chrono::high_resolution_clock::now() - time_now;
                    if(time_elapsed > update_interval){
                        time_now = std::chrono::high_resolution_clock::now();
                        umbrella_update_mesh(delta,3);
                        std::cout << " Updating Mesh " << std::endl;
                        // my_mesh.print_ds();
                    }
                }
                r.clear(glm::vec4(1.0, 1.0, 1.0, 1.0));

                camera.updateViewMatrix();

                Uint32 buttonState = SDL_GetMouseState(&xPos, &yPos);
                if( buttonState & SDL_BUTTON(SDL_BUTTON_LEFT) ) {
                    glm::vec4 pivot = glm::vec4(camera.lookAt.x, camera.lookAt.y, camera.lookAt.z, 1.0f);
                    glm::vec4 position = glm::vec4(camera.position.x, camera.position.y, camera.position.z, 1.0f);

                    float xAngle = (float)(lastxPos - xPos) * deltaAngleX;
                    float yAngle = (float)(lastyPos - yPos) * deltaAngleY;

                    float cosAngle = dot(camera.getViewDir(), camera.up);

                    if(cosAngle * signbit(deltaAngleY) > 0.99f)
                        deltaAngleY = 0.0f;

                    glm::mat4 rotationMatX(1.0f);
                    rotationMatX = glm::rotate(rotationMatX, xAngle, camera.up);
                    position = (rotationMatX * (position - pivot)) + pivot;

                    glm::mat4 rotationMatY(1.0f);
                    rotationMatY = glm::rotate(rotationMatY, yAngle, camera.getRightVector());
                    glm::vec3 finalPosition = (rotationMatY * (position - pivot)) + pivot;
                    camera.position = finalPosition;
                    camera.updateViewMatrix();
                }

                buttonState = SDL_GetMouseState(&xPos, &yPos);
                if( buttonState & SDL_BUTTON(SDL_BUTTON_RIGHT)) {
                    // Update camera parameters

                    float deltaY =  (float)(lastyPos - yPos) * 0.01f;
                    glm::mat4 dollyTransform = glm::mat4(1.0f);
                    dollyTransform = glm::translate(dollyTransform, normalize(camera.lookAt - camera.position) * deltaY);
                    glm::vec3 newCameraPosition = dollyTransform * glm::vec4(camera.position, 1.0f);
                    float newCameraFov = 2 * glm::atan(600.0f / (2 * deltaY)); // TODO Ask
                    
                    if(signbit(newCameraPosition.z) == signbit(camera.position.z)) {
                        camera.position = newCameraPosition;
                        camera.fov = newCameraFov; // TODO Ask
                        }
                }

                lastxPos = xPos;
                lastyPos = yPos;

                view = camera.getViewMatrix();
                
                r.setUniform(program, "model", model);
                r.setUniform(program, "view", view);
                r.setUniform(program, "projection", projection);
                r.setUniform(program, "lightPos", camera.position);
                r.setUniform(program, "viewPos", camera.position);
                r.setUniform(program, "lightColor", glm::vec3(1.0f, 1.0f, 1.0f));

                // r.setupFilledFaces();
                glm::vec3 red(1.0f, 0.0f, 0.0f);
                glm::vec3 green(0.0f, 1.0f, 0.0f);
                glm::vec3 blue(0.0f, 0.0f, 1.0f);
                glm::vec3 white(1.0f, 1.0f, 1.0f);
                r.setUniform(program, "ambientColor", 0.25f*white);
                r.setUniform(program, "intdiffuseColor", 0.75f*red);
                r.setUniform(program, "extdiffuseColor", 0.75f*blue);
                r.setUniform(program, "specularColor", 0.8f*white);
                r.setUniform(program, "phongExponent", 100.f);
                r.drawObject(object);

                r.setupWireFrame();
                glm::vec3 black(0.0f, 0.0f, 0.0f);
                r.setUniform(program, "ambientColor", black);
                r.setUniform(program, "intdiffuseColor", black);
                r.setUniform(program, "extdiffuseColor", black);
                r.setUniform(program, "specularColor", black);
                r.setUniform(program, "phongExponent", 0.f);
                r.drawEdges(wireframe);
                r.show();
            }
        }

    }
}
