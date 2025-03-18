#include <glm/gtc/matrix_transform.hpp>
#include <iostream>
#include <chrono>
#include <random>
#include <fstream>
#include <glad/gl.h>
#include <glm/glm.hpp>
#include <SDL2/SDL.h>
#include <string>

float a = 0.2;
float b = 0.2;
float c = 0.1;


// The parametric curve for the vase 
float radius_gen(float z) {
    return a + 
           b * glm::cos(glm::pi<float>() * (z - 0.2f)) +  // Shift bulge up
           c * glm::pow(z, 2) + 
           0.07f * glm::sin(2.0f * glm::pi<float>() * z) + // Gentle waviness
           0.06f * glm::exp(1.2f * z);  // Broaden the mouth
}



int main() {
    // We will use this file to generate an OBJ file

    float radius = 10.0f;
    float height = 20.0f;
    int slices = (int) (0.6f * height);

    // Base is xz plane
    // For a circle approxiamtion we will use a polygon with 36 sider

    // Create vertices
    std::vector<glm::vec3> vertices;
    // Create normals
    std::vector<glm::vec3> normals;

    // Add layer by layer
    int i = 0;
    while(i <= slices){
        // Use sin theta cos theta for x and z
        int j = 0;
        radius = radius_gen((float)i/(float)slices) * 15.0f;
        while(j < 36){
            float theta = 2.0f * 3.14f * (float)j / 36.0f;
            float x = radius * glm::cos(theta);
            float y = (float)i * height / (float)slices;
            float z = radius * glm::sin(theta);
            vertices.push_back(glm::vec3(x, y, z));
            normals.push_back(glm::normalize(glm::vec3(x, 0.0f, z)));
            j++;
        }
        i++;
    }

    // Now create the faces 
    std::vector<std::vector<int>> faces;
    // Faces will be created layer by layer in the form of rectangles
    i = 0;
    while(i < slices){
        int j = 0;
        while(j < 36){
            // Create a rectangle
            std::vector<int> rect;
            rect.push_back(i*36 + j);
            rect.push_back(i*36 + (j+1)%36);
            rect.push_back((i+1)*36 + (j+1)%36);
            rect.push_back((i+1)*36 + j);
            faces.push_back(rect);
            j++;
        }
        i++;
    }

    // Add bottom vertex
    vertices.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
    normals.push_back(glm::vec3(0.0f, -1.0f, 0.0f));
    int bottom_vertex = vertices.size() - 1;

    // Add the bottom face
    i = 0;
    while(i < 36){
        std::vector<int> bottom_face;
        bottom_face.push_back(i);
        bottom_face.push_back((i+1)%36);
        bottom_face.push_back(bottom_vertex);
        faces.push_back(bottom_face);
        i++;
    }

    // Add top vertex
    vertices.push_back(glm::vec3(0.0f, height, 0.0f));
    normals.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
    int top_vertex = vertices.size() - 1;
    // Add the top face
    i = 0;
    while(i < 36){
        std::vector<int> top_face;
        top_face.push_back((slices)*36 + i);
        top_face.push_back((slices)*36 + (i+1)%36);
        top_face.push_back(top_vertex);
        faces.push_back(top_face);
        i++;
    }

    // Now write to file ../../meshes/rajat.obj
    std::ofstream file("/Users/darkelixir/Desktop/COL781/A2/my_v/meshes/part_8.obj");
    for(const auto& vertex : vertices){
        file << "v " << vertex.x << " " << vertex.y << " " << vertex.z << "\n";
    }
    for(const auto& normal : normals){
        file << "vn " << normal.x << " " << normal.y << " " << normal.z << "\n";
    }
    for(const auto& face : faces){
        file << "f ";
        for(const auto& vertex : face){
            file << vertex + 1 << " ";
        }
        file << "\n";
    }
    file.close();

}