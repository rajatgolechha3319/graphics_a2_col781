#ifndef VIEWER_HPP
#define VIEWER_HPP

#include "../hw/hw.hpp"
#include "../../a2/a2.hpp"
namespace COL781 {
    namespace Viewer {

        class Camera {
        public:
            glm::vec3 position;
            glm::vec3 front;
            glm::vec3 up; 
            glm::vec3 lookAt;
            glm::mat4 viewMatrix;

            float cameraSpeed, yaw, pitch, lastX, lastY, fov, aspect;
            bool firstMouse;
            void initialize(float aspect);
            glm::mat4 getViewMatrix();
            glm::mat4 getProjectionMatrix();
            glm::vec3 getViewDir();
            glm::vec3 getRightVector();

            void setCameraView(glm::vec3 position_vector, glm::vec3 lookat_vector, glm::vec3 up_vector);
            void updateViewMatrix();
        };

        class Viewer {
        public:
            bool initialize(const std::string &title, int width, int height);
            void setMesh(int nv, int nt, int ne, const glm::vec3* vertices, const glm::ivec3* triangles, const glm::ivec2* edges, const glm::vec3* normals = nullptr);
            void setMesh_testing(int nv, int nt, int ne, const glm::vec3* vertices, const glm::ivec3* triangles, const glm::ivec2* edges, const glm::vec3* normals = nullptr);
            void catmull_clark();
            void setMesh_new(int nv, const glm::vec3* vertices, const std::vector<std::vector<int>> &poly_faces, const glm::vec3* normals = nullptr);
            void create_unit_rectangle(int m, int n);
            void create_sphere(int slices, int stacks);
            void create_cube(int m, int n, int o);
            void create_cube_new(int m, int n, int o);
            void create_noisy_cube(int m, int n, int o);
            void create_noisy_cube_new(int m, int n, int o);
            void umbrella_update_mesh(float delta, int iters);
            void load_obj_file(const std::string &filepath);
            void extrude(int face_idx, float d);
            void extrude_region(std::vector<int> face_idx, float d);
            void flip_normals();
            int get_closest_face(glm::vec3 p);
            void view(bool flag=false);
            private:
            COL781::OpenGL::Rasterizer r;
            mesh my_mesh;
            COL781::OpenGL::ShaderProgram program;
            COL781::OpenGL::Object object;
            COL781::OpenGL::Object wireframe;
            glm::mat4 stagetransform;
            Camera camera;
        };

    }
}

#endif
