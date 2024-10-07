/**
 * @file AI_functions.hpp
 * @author Hans Straile (straile@uni-potsdam.de)
 * @brief API for the AI/Machine Learning based chemistry surrogate model with functions to initialize a neural network and use it for training and inference via Keras for Python . 
 * @version 0.1
 * @date 01 Nov 2024
 *
 * This file implements the creation of a DHT by using the MPI
 * one-sided-communication. There is also the possibility to write or read data
 * from or to the DHT. In addition, the current state of the DHT can be written
 * to a file and read in again later.
 */

#ifndef AI_FUNCTIONS_H
#define AI_FUNCTIONS_H

#include <Python.h>
#include <string>

namespace poet {

int Python_Keras_setup(std::string functions_file_path);

int Python_Keras_load_model(std::string model_file_path);

PyObject* vector_to_numpy_array(const std::vector<std::vector<double>>& field);

std::vector<std::vector<double>> numpy_array_to_vector(PyObject* py_matrix);

std::vector<std::vector<double>> Python_keras_predict(std::vector<std::vector<double>> x, int batch_size);  

//int Python_keras_train(Field &x, Field &y, Field &x_val, Field &y_val, int batch_size, std::string pid)

} // namespace poet
#endif // AI_FUNCTIONS_HPP