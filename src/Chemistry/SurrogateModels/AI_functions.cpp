#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <Python.h>
#include <numpy/arrayobject.h>
#include <Eigen/Dense>
#include <thread>
#include <mutex>
#include <condition_variable>
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
  // Import numpy functions
  _import_array();
  FILE* fp = fopen(functions_file_path.c_str(), "r");
  int py_functions_initialized = PyRun_SimpleFile(fp, functions_file_path.c_str());
  fclose(fp);
  if (py_functions_initialized != 0) {
    PyErr_Print();
    throw std::runtime_error(std::string("AI surrogate Python functions could not be loaded." ) + 
                             "Are tensorflow and numpy installed?");
  }
  return py_functions_initialized;
}

/**
 * @brief Loads the user-supplied Keras model  
 * @param model_file_path Path to a .keras file that the user must supply as
 * a variable "model_file_path" in the R input script
 * @return 0 if function was succesful
 */
int Python_Keras_load_model(std::string model_file_path, std::string cuda_src_dir) {
  // Acquire the Python GIL
  PyGILState_STATE gstate = PyGILState_Ensure();

  // Initialize Keras model
  int py_model_loaded = PyRun_SimpleString(("model = initiate_model(\"" + 
                                            model_file_path + "\", \"" + 
                                            cuda_src_dir + "\")").c_str());
  if (py_model_loaded != 0) {
    PyErr_Print(); // Ensure that python errors make it to stdout
    throw std::runtime_error("Keras model could not be loaded from: " + model_file_path);
  }
  // Release the Python GIL
  PyGILState_Release(gstate);

  return py_model_loaded;
}

/**
 * @brief Converts the std::vector 2D matrix representation of a POET Field object to a numpy array
 * for use in the Python AI surrogate functions
 * @param field 2D-Matrix with the content of a Field object
 * @return Numpy representation of the input vector
 */
PyObject* vector_to_numpy_array(const std::vector<std::vector<double>>& field) {
  npy_intp dims[2] = {static_cast<npy_intp>(field[0].size()), 
                      static_cast<npy_intp>(field.size())};
  
  PyObject* np_array = PyArray_SimpleNew(2, dims, NPY_FLOAT64);
  double* data = static_cast<double*>(PyArray_DATA((PyArrayObject*)np_array));
  // write field data to numpy array
  for (size_t i = 0; i < field.size(); ++i) {
    for (size_t j = 0; j < field[i].size(); ++j) {
      data[j * field.size() + i] = field[i][j];
    }
  }
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
    std::cerr << "The model's output is not a numpy array." << std::endl;
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
std::vector<double> Python_Keras_predict(std::vector<std::vector<double>> x, int batch_size) {
  // Acquire the Python GIL
  PyGILState_STATE gstate = PyGILState_Ensure();
  
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
  Py_XDECREF(py_df_x);
  Py_XDECREF(args);
  Py_XDECREF(py_predictions);
  
  // Release the Python GIL
  PyGILState_Release(gstate);
  
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
Eigen::MatrixXd eigen_inference_batched(const Eigen::Ref<const Eigen::MatrixXd>& input_batch, const EigenModel& model) {
    Eigen::MatrixXd current_layer = input_batch;

    // Process all hidden layers
    for (size_t layer = 0; layer < model.weight_matrices.size() - 1; ++layer) {
        current_layer = (model.weight_matrices[layer] * current_layer);
        current_layer = current_layer.colwise() + model.biases[layer];
        current_layer = current_layer.array().max(0.0);
    }

    // Process output layer (without ReLU)
    size_t output_layer = model.weight_matrices.size() - 1;
    return (model.weight_matrices[output_layer] * current_layer).colwise() + model.biases[output_layer];
}

/**
 * @brief Uses the Eigen representation of the Keras model weights   for fast inference
 * @param x 2D-Matrix with the content of a Field object
 * @param batch_size size for mini-batches that are used in the Keras model.predict() method
 * @return Predictions that the neural network made from the input values x. The predictions are
 * represented as a vector similar to the representation from the Field.AsVector() method
 */
std::vector<double> Eigen_predict(const EigenModel& model, std::vector<std::vector<double>> x, int batch_size,
                                  std::mutex* Eigen_model_mutex) {
  // Convert input data to Eigen matrix
  const int num_samples = x[0].size();
  const int num_features = x.size();
  Eigen::MatrixXd full_input_matrix(num_features, num_samples);

  for (int i = 0; i < num_samples; ++i) {
    for (int j = 0; j < num_features; ++j) {
      full_input_matrix(j, i) = x[j][i];
    }
  }

  std::vector<double> result;
  result.reserve(num_samples * num_features);
  
  if (num_features != model.weight_matrices[0].cols()) {
    throw std::runtime_error("Input data size " + std::to_string(num_features) + \
      " does not match model input layer of size " + std::to_string(model.weight_matrices[0].cols()));
  }
  int num_batches = std::ceil(static_cast<double>(num_samples) / batch_size);

  Eigen_model_mutex->lock();
  for (int batch = 0; batch < num_batches; ++batch) {
    int start_idx = batch * batch_size;
    int end_idx = std::min((batch + 1) * batch_size, num_samples);
    int current_batch_size = end_idx - start_idx;
    // Extract the current input data batch
    Eigen::MatrixXd batch_data(num_features, current_batch_size);
    batch_data = full_input_matrix.block(0, start_idx, num_features, current_batch_size);
    // Predict
    batch_data = eigen_inference_batched(batch_data, model);

    result.insert(result.end(), batch_data.data(), batch_data.data() + batch_data.size());
  }
  Eigen_model_mutex->unlock();
  return result;
}


/**
 * @brief Appends data from one matrix (column major std::vector<std::vector<double>>) to another
 * @param training_data_buffer Matrix that the values are appended to
 * @param new_values Matrix that is appended
 */
void training_data_buffer_append(std::vector<std::vector<double>>& training_data_buffer, 
                                 std::vector<std::vector<double>>& new_values) {
  // Initialize training data buffer if empty
  if (training_data_buffer.size() == 0) {
    training_data_buffer = new_values;
  } else { // otherwise append
    for (int col = 0; col < training_data_buffer.size(); col++) {
      training_data_buffer[col].insert(training_data_buffer[col].end(),
                                        new_values[col].begin(), new_values[col].end());
    }
  }
}


/**
 * @brief Uses the Python environment with Keras' default functions to train the model
 * @param x Training data features
 * @param y Training data targets
 * @param params Global runtime paramters
 */
void Python_keras_train(std::vector<std::vector<double>> x, std::vector<std::vector<double>> y,
                        const RuntimeParameters& params) {
  // Prepare data for python
  PyObject* py_df_x = vector_to_numpy_array(x);
  PyObject* py_df_y = vector_to_numpy_array(y);

  // Get the model and training function from the global python interpreter
  PyObject* py_main_module = PyImport_AddModule("__main__");
  PyObject* py_global_dict = PyModule_GetDict(py_main_module);
  PyObject* py_keras_model = PyDict_GetItemString(py_global_dict, "model");
  PyObject* py_training_function = PyDict_GetItemString(py_global_dict, "training_step");
  
  // Build the function arguments as four python objects and an integer
  PyObject* args = Py_BuildValue("(OOOiis)", 
      py_keras_model, py_df_x, py_df_y, params.batch_size, params.training_epochs,
      params.save_model_path.c_str());
  

  // Call the Python training function
  PyObject *py_rv = PyObject_CallObject(py_training_function, args);

  // Clean up
  PyErr_Print(); 
  Py_DECREF(py_df_x);
  Py_DECREF(py_df_y);
  Py_DECREF(args);
}

/**
 * @brief Function for threadsafe parallel training and weight updating.
 * The function waits conditionally until the training data buffer is full.
 * It then clears the buffer and starts training, after training it writes the new weights to 
 * the Eigen model.  
 * @param Eigen_model Pointer to the EigenModel struct that will be updates with new weights
 * @param Eigen_model_mutex Mutex to ensure threadsafe access to the EigenModel struct
 * @param training_data_buffer Pointer to the Training data struct with which the model is trained
 * @param training_data_buffer_mutex Mutex to ensure threadsafe access to the training data struct
 * @param training_data_buffer_full Conditional waiting variable with wich the main thread signals
 * when a training run can start
 * @param start_training Conditional waiting predicate to mitigate against spurious wakeups
 * @param end_training Signals end of program to wind down thread gracefully 
 * @param params Global runtime paramters
 * @return 0 if function was succesful
 */
void parallel_training(EigenModel* Eigen_model,
                       std::mutex* Eigen_model_mutex,
                       TrainingData* training_data_buffer,
                       std::mutex* training_data_buffer_mutex,
                       std::condition_variable* training_data_buffer_full,
                       bool* start_training, bool* end_training,
                       const RuntimeParameters& params) {
  while (true) {
    // Conditional waiting:
    // - Sleeps until a signal arrives on training_data_buffer_full
    // - Releases the lock on training_data_buffer_mutex while sleeping
    // - Lambda function with start_training checks if it was a spurious wakeup
    // - Reaquires the lock on training_data_buffer_mutex after waking up
    // - If start_training has been set to true while the thread was active, it does NOT
    //   wait for a signal on training_data_buffer_full but starts the next round immediately.
    std::unique_lock<std::mutex> lock(*training_data_buffer_mutex);
    training_data_buffer_full->wait(lock, [start_training] { return *start_training;});
    
    //hier nochmal training_data_buffer_mutex lock/lock test?
    
    // Return if program is about to end
    if (*end_training) {
      return;
    }
    
    // Reset the  waiting predicate
    *start_training = false;
    // Get the necessary training data
    std::cout << "AI: Training thread: Getting training data" << std::endl;
    // Initialize training data input and targets
    std::vector<std::vector<double>> inputs(9);
    std::vector<std::vector<double>> targets(9);
    for (int col = 0; col < training_data_buffer->x.size(); col++) {
      // Copy data from the front of the buffer to the training inputs
      std::copy(training_data_buffer->x[col].begin(),
                training_data_buffer->x[col].begin() + params.training_data_size, 
                std::back_inserter(inputs[col]));
      // Remove copied data from the front of the buffer
      training_data_buffer->x[col].erase(training_data_buffer->x[col].begin(),
                                         training_data_buffer->x[col].begin() + params.training_data_size);

      // Copy data from the front of the buffer to the training targets
      std::copy(training_data_buffer->y[col].begin(),
                training_data_buffer->y[col].begin() + params.training_data_size, 
                std::back_inserter(targets[col]));
      // Remove copied data from the front of the buffer
      training_data_buffer->y[col].erase(training_data_buffer->y[col].begin(),
                                         training_data_buffer->y[col].begin() + params.training_data_size);
    }
    //update number of training runs
    training_data_buffer->n_training_runs += 1;
    // Unlock the training_data_buffer_mutex 
    lock.unlock();

    std::cout << "AI: Training thread: Start training" << std::endl;

    // Acquire the Python GIL
    PyGILState_STATE gstate = PyGILState_Ensure();

    // Start training
    Python_keras_train(inputs, targets, params);
    
    if (!params.use_Keras_predictions) {
      std::cout << "AI: Training thread: Update shared model weights" << std::endl;
      std::vector<std::vector<std::vector<double>>> cpp_weights = Python_Keras_get_weights();
      Eigen_model_mutex->lock();
      update_weights(Eigen_model, cpp_weights);
      Eigen_model_mutex->unlock();
    }
    
    // Release the Python GIL
    PyGILState_Release(gstate);
    std::cout << "AI: Training thread: Finished training, waiting for new data" << std::endl;
  } 
}

std::thread python_train_thread;
/**
 * @brief Starts a thread for parallel training and weight updating. This Wrapper function
 * ensures, that the main POET program can be built without pthread support if the AI
 * surrogate functions are disabled during compilation. 
 * @param Eigen_model Pointer to the EigenModel struct that will be updates with new weights
 * @param Eigen_model_mutex Mutex to ensure threadsafe access to the EigenModel struct
 * @param training_data_buffer Pointer to the Training data struct with which the model is trained
 * @param training_data_buffer_mutex Mutex to ensure threadsafe access to the training data struct
 * @param training_data_buffer_full Conditional waiting variable with wich the main thread signals
 * when a training run can start
 * @param start_training Conditional waiting predicate to mitigate against spurious wakeups
 * @param end_training Signals end of program to wind down thread gracefully 
 * @param params Global runtime paramters
 * @return 0 if function was succesful
 */
int Python_Keras_training_thread(EigenModel* Eigen_model,
                                 std::mutex* Eigen_model_mutex,
                                 TrainingData* training_data_buffer,
                                 std::mutex* training_data_buffer_mutex,
                                 std::condition_variable* training_data_buffer_full,
                                 bool* start_training, bool* end_training,
                                 const RuntimeParameters& params) {

  PyThreadState *_save = PyEval_SaveThread();
  python_train_thread = std::thread(parallel_training, Eigen_model, Eigen_model_mutex,
                                    training_data_buffer, training_data_buffer_mutex,
                                    training_data_buffer_full, start_training, end_training,
                                    params);
  return 0;
}


/**
 * @brief Updates the EigenModels weigths and biases from the weight vector
 * @param model Pounter to an EigenModel struct
 * @param weights Cector of model weights from keras as returned by Python_Keras_get_weights()
 */
void update_weights(EigenModel* model,
                       const std::vector<std::vector<std::vector<double>>>& weights) {
  size_t num_layers = weights.size() / 2;
  for (size_t i = 0; i < weights.size(); i += 2) {
      // Fill current weight matrix
      size_t rows = weights[i][0].size();
      size_t cols = weights[i].size();
      for (size_t j = 0; j < cols; ++j) {
          for (size_t k = 0; k < rows; ++k) {
              model->weight_matrices[i / 2](k, j) = weights[i][j][k];
          }
      }
      // Fill bias vector
      size_t bias_size = weights[i + 1][0].size();
      for (size_t j = 0; j < bias_size; ++j) {
          model->biases[i / 2](j) = weights[i + 1][0][j];
      }
  }
}


/**
 * @brief Converts the weights and biases from the Python Keras model to C++ vectors
 * @return A vector containing the model weights and biases 
 */
std::vector<std::vector<std::vector<double>>> Python_Keras_get_weights() {
  // Acquire the Python GIL
  PyGILState_STATE gstate = PyGILState_Ensure();

  PyObject* py_main_module = PyImport_AddModule("__main__");
  PyObject* py_global_dict = PyModule_GetDict(py_main_module);
  PyObject* py_keras_model = PyDict_GetItemString(py_global_dict, "model");
  PyObject* py_get_weights_function = PyDict_GetItemString(py_global_dict, "get_weights");
  PyObject* args = Py_BuildValue("(O)", py_keras_model);
  
  // Call Python function
  PyObject* py_weights_list = PyObject_CallObject(py_get_weights_function, args);
  
  if (!py_weights_list) {
      PyErr_Print();
      throw std::runtime_error("Failed to get weights from Keras model");
  }

  // Container for the extracted weights
  std::vector<std::vector<std::vector<double>>> cpp_weights;

  // Iterate through the layers (weights and biases)
  Py_ssize_t num_layers = PyList_Size(py_weights_list);
  for (Py_ssize_t i = 0; i < num_layers; ++i) {
    PyObject* py_weight_array = PyList_GetItem(py_weights_list, i);
    if (!PyArray_Check(py_weight_array)) {
      throw std::runtime_error("Weight is not a NumPy array.");
    }

    PyArrayObject* weight_np = reinterpret_cast<PyArrayObject*>(py_weight_array);
    int dtype = PyArray_TYPE(weight_np);

    if (PyArray_NDIM(weight_np) == 2) {
      int num_rows = PyArray_DIM(weight_np, 0);
      int num_cols = PyArray_DIM(weight_np, 1);
      
      std::vector<std::vector<double>> weight_matrix(num_rows, std::vector<double>(num_cols));
      // Handle different precision settings
      if (dtype == NPY_FLOAT32) {
        float* weight_data_float = static_cast<float*>(PyArray_DATA(weight_np));
        for (int r = 0; r < num_rows; ++r) {
          for (int c = 0; c < num_cols; ++c) {
            weight_matrix[r][c] = static_cast<double>(weight_data_float[r * num_cols + c]);
          }
        }
      } else if (dtype == NPY_DOUBLE) {
        double* weight_data_double = static_cast<double*>(PyArray_DATA(weight_np));
        for (int r = 0; r < num_rows; ++r) {
          for (int c = 0; c < num_cols; ++c) {
            weight_matrix[r][c] = weight_data_double[r * num_cols + c];
          }
        }
      } else {
          throw std::runtime_error("Unsupported data type for weights. Must be NPY_FLOAT32 or NPY_DOUBLE.");
      }
      cpp_weights.push_back(weight_matrix);

    } else if (PyArray_NDIM(weight_np) == 1) {
      // 1D bias vectors
      int num_elements = PyArray_DIM(weight_np, 0);
      std::vector<std::vector<double>> bias_vector(1, std::vector<double>(num_elements));

      // Handle different precision settings
      if (dtype == NPY_FLOAT32) {
        float* bias_data_float = static_cast<float*>(PyArray_DATA(weight_np));
        for (int j = 0; j < num_elements; ++j) {
          bias_vector[0][j] = static_cast<double>(bias_data_float[j]);
        }
      } else if (dtype == NPY_DOUBLE) {
        double* bias_data_double = static_cast<double*>(PyArray_DATA(weight_np));
        for (int j = 0; j < num_elements; ++j) {
          bias_vector[0][j] = bias_data_double[j];
        }
      } else {
        throw std::runtime_error("Unsupported data type for biases. Must be NPY_FLOAT32 or NPY_DOUBLE.");
      }
      cpp_weights.push_back(bias_vector);
    } else {
      throw std::runtime_error("Unsupported weight dimension. Only 1D and 2D arrays are supported.");
    }
  }
  // Clean up
  Py_DECREF(py_weights_list);
  Py_DECREF(args);
  // Release Python GIL
  PyGILState_Release(gstate);

  return cpp_weights;
}


/**
 * @brief Joins the training thread and winds down the Python environment gracefully
 * @param Eigen_model_mutex Mutex to ensure threadsafe access to the EigenModel struct
 * @param training_data_buffer_mutex Mutex to ensure threadsafe access to the training data struct
 * @param training_data_buffer_full Conditional waiting variable with wich the main thread signals
 * when a training run can start
 * @param start_training Conditional waiting predicate to mitigate against spurious wakeups
 * @param end_training Signals end of program to wind down thread gracefully */
void Python_finalize(std::mutex* Eigen_model_mutex, std::mutex* training_data_buffer_mutex,
                     std::condition_variable* training_data_buffer_full,
                     bool* start_training, bool* end_training) {
  // Define training as over
  *end_training = true;
  // Wake up and join training thread
  *start_training = true;
  training_data_buffer_full->notify_one();
  
  if (python_train_thread.joinable()) {
    python_train_thread.join();
  }

  // Acquire the Python GIL
  PyGILState_STATE gstate = PyGILState_Ensure();
  //Finalize Python
  Py_FinalizeEx();
}

} //namespace poet