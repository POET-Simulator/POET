#include "AI_functions.hpp"
#include <iostream>
#include <string>
#include <cstring>
#include <Python.h>

using namespace std;

namespace poet {

void Python_Keras_setup(std::string SRC_DIR, std::string model_file_path) {
    Py_Initialize();  // Initialize the Python interpreter
    std::cout << SRC_DIR + 
                "/src/Chemistry/SurrogateModels/AI_Python_functions/keras_AI_surrogate.py" 
              << std::endl;
    std::string python_keras_file = SRC_DIR + "/src/Chemistry/SurrogateModels/AI_Python_functions/keras_AI_surrogate.py"; 
    FILE* fp = fopen(python_keras_file.c_str(), "r");
    PyRun_SimpleFile(fp, python_keras_file.c_str());
    PyRun_SimpleString(("model = initiate_model(\"" + model_file_path + "\")").c_str());
    PyErr_Print(); // Ensure that python errors make it to stdout
    fclose(fp);
}

}