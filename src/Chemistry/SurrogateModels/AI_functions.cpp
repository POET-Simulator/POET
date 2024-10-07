#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <Python.h>
#include <numpy/arrayobject.h>
#include <Eigen/Dense>
#include "AI_functions.hpp"


using namespace std;

namespace poet {


/**
 * @brief Loads the Python interpreter and functions  
 * @param functions_file_path Path to the Python file where the AI surrogate
 * functions are defined
 * @return 0 if function was succesful
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
 * @return 0 if function was succesful
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
 * @return Numpy representation of the input vector
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
 * @brief Converts a Pyton matrix object to a std::vector vector
 * @param py_matrix Pyobject that must be a 2D matrix
 * @result Vector that can be used similar to the return value of the Field object
 * Field.AsVector() method.
 */
std::vector<double> numpy_array_to_vector(PyObject* py_array) {
  std::vector<double> result;
  if (!PyArray_Check(py_array)) {
    std::cerr << "The model's predictions are not a numpy array." << std::endl;
    return result;
  }
  // Cast generic PyObject to PyArrayObject
  PyArrayObject* np_array = reinterpret_cast<PyArrayObject*>(py_array);
  
  // Get shape
  int numDims = PyArray_NDIM(np_array);
  npy_intp* shape = PyArray_SHAPE(np_array);
  if (numDims != 2) {
    std::cerr << "The model's predictions are not a 2D matrix." << std::endl;
    return result;
  }

  // Copy data into std::vector format
  double* data = static_cast<double*>(PyArray_DATA(np_array));
  npy_intp size = PyArray_SIZE(np_array);
  result.resize(size);
  std::copy(data, data + size, result.begin());
  return result;
}

/**
 * @brief Uses the Python Keras functions to calculate predictions from a neural network.
 * @param x 2D-Matrix with the content of a Field object
 * @param batch_size size for mini-batches that are used in the Keras model.predict() method
 * @return Predictions that the neural network made from the input values x. The predictions are
 * represented as a vector similar to the representation from the Field.AsVector() method
 */
std::vector<double> Python_keras_predict(std::vector<std::vector<double>> x, int batch_size) {  
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
  std::vector<double> predictions = numpy_array_to_vector(py_predictions); 
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

/**
 * @brief Uses Eigen for fast inference with the weights and biases of a neural network 
 * @param input_batch Batch of input data that must fit the size of the neural networks input layer
 * @param model Struct of aligned Eigen vectors that hold the neural networks weights and biases.
 * Only supports simple fully connected feed forward networks.
 * @return The batch of predictions made with the neural network weights and biases and the data
 * in input_batch
 */
Eigen::MatrixXd eigen_inference_batched(const Eigen::Ref<Eigen::MatrixXd>& input_batch, const EigenModel& model) {
    Eigen::MatrixXd current_layer = input_batch;
    std::cout << "input_batch ROWS " << current_layer.rows() << std::endl;
    std::cout << "input_batch COLS " << current_layer.cols() << std::endl;
    // Process all hidden layers
    for (size_t layer = 0; layer < model.weight_matrices.size() - 1; ++layer) {

        std::cout << "WEIGHTS LAYER " << layer << std::endl;
        std::cout << "WEIGHTS ROWS " << model.weight_matrices[layer].rows() << std::endl;
        std::cout << "WEIGHTS COLS " << model.weight_matrices[layer].cols() << std::endl;
        std::cout << "BIASES SIZE " << model.biases[layer].size() << std::endl;

        current_layer = (model.weight_matrices[layer] * current_layer);
        std::cout << "MULTIPLIED" << std::endl;
        current_layer = current_layer.colwise() + model.biases[layer];
        current_layer = current_layer.array().max(0.0);
    }

    // Process output layer (without ReLU)
    size_t output_layer = model.weight_matrices.size() - 1;
    return (model.weight_matrices[output_layer] * current_layer).colwise() + model.biases[output_layer];
}
/**
 * @brief Converts the weights and biases from a Python Keras model to Eigen matrices
 * @return The struct containing the model weights and biases as aligned Eigen matrices 
 */
#include <iostream>

EigenModel Python_Keras_get_weights_as_Eigen() {
    EigenModel eigen_model;

    PyObject* py_main_module = PyImport_AddModule("__main__");
    if (!py_main_module) throw std::runtime_error("Failed to import Python main module");

    PyObject* py_global_dict = PyModule_GetDict(py_main_module);
    if (!py_global_dict) throw std::runtime_error("Failed to get global dictionary");

    PyObject* py_keras_model = PyDict_GetItemString(py_global_dict, "model");
    if (!py_keras_model) throw std::runtime_error("Failed to get Keras model from Python dictionary");

    PyObject* py_get_weights_function = PyDict_GetItemString(py_global_dict, "get_weights");
    if (!py_get_weights_function) throw std::runtime_error("Failed to get Keras get_weights function");

    PyObject* args = Py_BuildValue("(O)", py_keras_model);
    if (!args) throw std::runtime_error("Failed to create argument tuple");

    // Call Python function
    PyObject* py_weights_list = PyObject_CallObject(py_get_weights_function, args);
    Py_DECREF(args); // Cleanup args

    if (!py_weights_list) {
        PyErr_Print(); // Print Python error
        throw std::runtime_error("Failed to get weights from Keras model");
    }

    Py_ssize_t num_layers = PyList_Size(py_weights_list);
    for (Py_ssize_t i = 0; i < num_layers; i += 2) {
        // Get weight matrix
        PyObject* weight_array = PyList_GetItem(py_weights_list, i);
        if (!weight_array) throw std::runtime_error("Failed to get weight array");
        if (!PyArray_Check(weight_array)) throw std::runtime_error("Weight array is not a NumPy array");

        PyArrayObject* weight_np = reinterpret_cast<PyArrayObject*>(weight_array);
        if (PyArray_NDIM(weight_np) != 2) throw std::runtime_error("Weight array is not 2-dimensional");

        // Check data type
        int dtype = PyArray_TYPE(weight_np);
        if (dtype == NPY_FLOAT32) {
          std::cout << "DTYPE: SINGLE" << std::endl;
        }
        if (dtype == NPY_DOUBLE) {
          std::cout << "DTYPE: DOUBLE" << std::endl;
        }
        if (dtype != NPY_FLOAT32 && dtype != NPY_DOUBLE) {
            throw std::runtime_error("Unsupported data type for weight array");
        }

        // Cast the data correctly depending on the type
        void* weight_data = PyArray_DATA(weight_np);
        if (!weight_data) throw std::runtime_error("Failed to get weight data pointer");

        int num_rows = PyArray_DIM(weight_np, 0);  // NumPy's row count
        int num_cols = PyArray_DIM(weight_np, 1);  // NumPy's column count

        Eigen::MatrixXd weight_matrix; // Use double-precision matrix

        if (dtype == NPY_FLOAT32) {
            // Handle float32 array from NumPy
            float* weight_data_float = static_cast<float*>(weight_data);
            weight_matrix = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
                weight_data_float, num_rows, num_cols).cast<double>().transpose();
        } else if (dtype == NPY_DOUBLE) {
            // Handle double (float64) array from NumPy
            double* weight_data_double = static_cast<double*>(weight_data);
            weight_matrix = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
                weight_data_double, num_rows, num_cols).transpose();
        }

        // Get bias vector
        PyObject* bias_array = PyList_GetItem(py_weights_list, i + 1);
        if (!bias_array) throw std::runtime_error("Failed to get bias array");
        if (!PyArray_Check(bias_array)) throw std::runtime_error("Bias array is not a NumPy array");

        PyArrayObject* bias_np = reinterpret_cast<PyArrayObject*>(bias_array);
        if (PyArray_NDIM(bias_np) != 1) throw std::runtime_error("Bias array is not 1-dimensional");

        // Check bias data type and cast accordingly
        int bias_dtype = PyArray_TYPE(bias_np);
        void* bias_data = PyArray_DATA(bias_np);
        if (!bias_data) throw std::runtime_error("Failed to get bias data pointer");

        int bias_size = PyArray_DIM(bias_np, 0);
        Eigen::VectorXd bias_vector; // Use double-precision vector

        if (bias_dtype == NPY_FLOAT32) {
            float* bias_data_float = static_cast<float*>(bias_data);
            bias_vector = Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, 1>>(bias_data_float, bias_size).cast<double>();
        } else if (bias_dtype == NPY_DOUBLE) {
            double* bias_data_double = static_cast<double*>(bias_data);
            bias_vector = Eigen::Map<Eigen::VectorXd>(bias_data_double, bias_size);
        }

        // Add to EigenModel
        eigen_model.weight_matrices.push_back(weight_matrix);
        eigen_model.biases.push_back(bias_vector);
    }

    // Clean up
    Py_DECREF(py_weights_list); // Clean up Python list object
    Py_DECREF(py_main_module); // Clean up Python module object

    return eigen_model;
}





/**
 * @brief Gets the model weights and biases of the Keras model and uses Eigen for fast inference
 * @param x 2D-Matrix with the content of a Field object
 * @param batch_size size for mini-batches that are used in the Keras model.predict() method
 * @return Predictions that the neural network made from the input values x. The predictions are
 * represented as a vector similar to the representation from the Field.AsVector() method
 */
std::vector<double> Eigen_predict(const EigenModel& model, std::vector<std::vector<double>> x, int batch_size) {
  // Convert input data to Eigen matrix
  const int num_samples = x[0].size();
  const int num_features = x.size();
  
  std::cout << "num_samples " << num_samples << std::endl;
  std::cout << "num_features " << num_features << std::endl;

  Eigen::MatrixXd full_input_matrix(num_features, num_samples);
  for (int i = 0; i < num_samples; ++i) {
    for (int j = 0; j < num_features; ++j) {
      full_input_matrix(j, i) = x[j][i];
      if (i < 6) {std::cout << full_input_matrix.coeff(j, i) << "  ";}
    }
    if (i < 6) {std::cout << std::endl;}
  }
    

  std::cout << "Eigen rows " << full_input_matrix.rows() << std::endl;
  std::cout << "Eigen cols " << full_input_matrix.cols() << std::endl;

    std::vector<double> result;
    result.reserve(num_samples * model.biases.back().size());
    
    int num_batches = std::ceil(static_cast<double>(num_samples) / batch_size);
    for (int batch = 0; batch < num_batches; ++batch) {
        int start_idx = batch * batch_size;
        int end_idx = std::min((batch + 1) * batch_size, num_samples);
        int current_batch_size = end_idx - start_idx;
        std::cout << "CREATE INPUT BATCH" << std::endl;
        // Extract the current batch
        Eigen::MatrixXd input_batch(num_features, current_batch_size); 
        input_batch = full_input_matrix.block(0, start_idx, num_features, current_batch_size);
        // Perform inference on the batch
        Eigen::MatrixXd output_batch(num_features, current_batch_size);
        output_batch = eigen_inference_batched(input_batch, model);
        for (int i = 0; i < num_samples; ++i) {
          for (int j = 0; j < num_features; ++j) {
            if (i < 6) {std::cout << output_batch.coeff(j, i) << "  ";}
          }
          if (i < 6) {std::cout << std::endl;}
        }
        
        //result.insert(result.end(), output_batch.data(), output_batch.data() + output_batch.size());
    }
    
    return result;
}
}