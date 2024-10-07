#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <Python.h>
#include <numpy/arrayobject.h>
#include "AI_functions.hpp"


using namespace std;

namespace poet {


/**
 * @brief Loads the Python interpreter and functions  
 * @param functions_file_path Path to the Python file where the AI surrogate
 * functions are defined
 */
int Python_Keras_setup(std::string functions_file_path) {
  // Initialize Python functions
  Py_Initialize();
  FILE* fp = fopen(functions_file_path.c_str(), "r");
  int py_functions_initialized = PyRun_SimpleFile(fp, functions_file_path.c_str());
  fclose(fp);
  if (py_functions_initialized != 0) {
    PyErr_Print();
    throw std::runtime_error("Python functions could not be initialized");
  }
  return py_functions_initialized;
}

/**
 * @brief Loads the user-supplied Keras model  
 * @param model_file_path Path to a .keras file that the user must supply as
 * a variable "model_file_path" in the R input script
 */
int Python_Keras_load_model(std::string model_file_path) {
  // Initialize Keras model
  int py_model_loaded = PyRun_SimpleString(("model = initiate_model(\"" + 
                                            model_file_path + "\")").c_str());
  if (py_model_loaded != 0) {
    PyErr_Print(); // Ensure that python errors make it to stdout
    throw std::runtime_error("Keras model could not be loaded from: " + model_file_path);
  }
  return py_model_loaded;
}

/**
 * @brief Converts the std::vector 2D matrix representation of a POET Field object to a numpy array
 * for use in the Python AI surrogate functions
 * @param field 2D-Matrix with the content of a Field object
 */
PyObject* vector_to_numpy_array(const std::vector<std::vector<double>>& field) {
  // import numpy and numpy array API
  PyObject* numpy = PyImport_ImportModule("numpy");
  import_array();
  npy_intp dims[2] = {static_cast<npy_intp>(field[0].size()), 
                      static_cast<npy_intp>(field.size())};
  PyObject* np_array = PyArray_SimpleNew(2, dims, NPY_FLOAT64);
  _Float64* data = static_cast<_Float64*>(PyArray_DATA((PyArrayObject*)np_array));
  // write field data to numpy array
  for (size_t i = 0; i < field.size(); ++i) {
    for (size_t j = 0; j < field[i].size(); ++j) {
      data[j * field.size() + i] = field[i][j];
    }
  }
  // Clean up
  Py_DECREF(numpy);
  return np_array;
}

/**
 * @brief Converts a Pyton matrix object to a std::vector vector 2D matrix  
 * @param py_matrix Pyobject that must be a 2D matrix
 */
std::vector<std::vector<double>> numpy_array_to_vector(PyObject* pyArray) {
  std::vector<std::vector<double>> result;
  if (!PyArray_Check(pyArray)) {
    std::cerr << "The model's predictions are not a numpy array." << std::endl;
    return result;
  }
  // Cast generic PyObject to PyArrayObject
  PyArrayObject* npArray = reinterpret_cast<PyArrayObject*>(pyArray);
  
  // Get shape
  int numDims = PyArray_NDIM(npArray);
  npy_intp* shape = PyArray_SHAPE(npArray);

  if (numDims != 2) {
    std::cerr << "The model's predictions are not a 2D matrix." << std::endl;
    return result;
  }

  // Copy data into std::vector format
  double* data = static_cast<double*>(PyArray_DATA(npArray));
  result.resize(shape[0]);
  for (npy_intp i = 0; i < shape[0]; ++i) {
    result[i].resize(shape[1]);
    for (npy_intp j = 0; j < shape[1]; ++j) {
      result[i][j] = data[i * shape[1] + j];
    }
  }
  return result;
}


/**
 * @brief Converts the 2D-Matrix represantation of a POET Field object to a numpy array
 * for use in the Python AI surrogate functions
 * @param field 2D-Matrix with the content of a Field object
 */
std::vector<std::vector<double>> Python_keras_predict(std::vector<std::vector<double>> x, int batch_size) {  
  // Prepare data for python
  PyObject* py_df_x = vector_to_numpy_array(x);
  
  // Get the model and training function from the global python interpreter
  PyObject* py_main_module = PyImport_AddModule("__main__");
  PyObject* py_global_dict = PyModule_GetDict(py_main_module);
  PyObject* py_keras_model = PyDict_GetItemString(py_global_dict, "model");
  PyObject* py_inference_function = PyDict_GetItemString(py_global_dict, "prediction_step");
   // Build the function arguments as four python objects and an integer
  PyObject* args = Py_BuildValue("(OOi)",
      py_keras_model, py_df_x, batch_size);
  
  // Call the Python training function
  PyObject *py_predictions = PyObject_CallObject(py_inference_function, args);

  // Check py_rv and return as 2D vector
  std::vector<std::vector<double>> predictions = numpy_array_to_vector(py_predictions); 
  // Clean up
  PyErr_Print(); // Ensure that python errors make it to stdout 
  Py_DECREF(py_df_x);
  Py_DECREF(py_main_module);
  Py_DECREF(py_global_dict);
  Py_DECREF(py_keras_model);
  Py_DECREF(py_inference_function);
  Py_DECREF(args);
  Py_DECREF(py_predictions);
  return predictions;
} 

}