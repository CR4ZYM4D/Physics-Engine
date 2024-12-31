#include <GLFW/glfw3.h>
#include <iostream>
#include "Included files/particles.h"

using namespace std;
int main(void)
{
    GLFWwindow* window;
    Vector3D v1 = Vector3D(3, 1, 1);
    real s = 9;
    v1 *= s;
    cout << v1.x << endl ;
    Vector3D v2 = Vector3D(3,1,1);
    cout << v1.X(v2).z;
    //Initialize the library 
    if (!glfwInit())
        return -1;

     //Create a windowed mode window and its OpenGL context 
    window = glfwCreateWindow(900, 700, "Hello World", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    // Make the window's context current 
    glfwMakeContextCurrent(window);

    // Loop until the user closes the window 
    while (!glfwWindowShouldClose(window))
    {
        // Render here 
        glClear(GL_COLOR_BUFFER_BIT);

        // Swap front and back buffers 
        glfwSwapBuffers(window);

        // Poll for and process events 
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}