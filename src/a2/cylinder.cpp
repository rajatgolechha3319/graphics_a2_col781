#include <glm/gtc/matrix_transform.hpp>
#include <iostream>
#include <chrono>
#include <random>
#include <glad/gl.h>
#include <glm/glm.hpp>
#include <SDL2/SDL.h>
#include <string>


int main() {
    // We will use this file to generate an OBJ file

    float radius = 1.0f;
    int slices = 20;
    float height = 20.0f;


    // Base is xz plane
    // For a circle approxiamtion we will use a polygon with 36 sider

    // Create vertices
    std::vector<glm::vec3> vertices;

    // Add layer by layer
    int i = 0;
    while(i <= slices){
        // Use sin theta cos theta for x and z
        int j = 0;
        while(j < 36){
            float theta = 2.0f * 3.14f * (float)j / 36.0f;
            float x = radius * glm::cos(theta);
            float z = radius * glm::sin(theta);
            float y = (float)i * height / (float)slices;
            vertices.push_back(glm::vec3(x, y, z));
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

    // Now write to file ../../meshes/rajat.obj
    std::ofstream file("../../meshes/rajat.obj");
    for(const auto& vertex : vertices){
        file << "v " << vertex.x << " " << vertex.y << " " << vertex.z << "\n";
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