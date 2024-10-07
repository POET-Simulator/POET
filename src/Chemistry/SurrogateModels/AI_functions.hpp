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

#include <string>

namespace poet {

void Python_Keras_setup(std::string SRC_DIR, std::string model_file_path);

} // namespace poet
#endif // AI_FUNCTIONS_HPP