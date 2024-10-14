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
#include <vector>
#include "poet.hpp"
#include "Base/RInsidePOET.hpp"
#include "Init/InitialList.hpp"

// PhreeqC definition of pi clashes with Eigen macros so we have to temporarily undef it 
#pragma push_macro("pi")
#undef pi
#include <Eigen/Dense> 
#pragma pop_macro("pi")

namespace poet {
// Define an aligned allocator for std::vector
template<typename T>
using aligned_vector = std::vector<T, Eigen::aligned_allocator<T>>;
// Define a structure to hold the weights in Eigen matrices
struct EigenModel {
  aligned_vector<Eigen::MatrixXd> weight_matrices;
  aligned_vector<Eigen::VectorXd> biases;
};

struct TrainingData {
  std::vector<std::vector<double>> x;
  std::vector<std::vector<double>> y;
};

// Ony declare the actual functions if flag is set 
#ifdef USE_AI_SURROGATE

int Python_Keras_setup(std::string functions_file_path);

void set_ai_surrogate_runtime_params(RInsidePOET& R, RuntimeParameters& params, InitialList& init_list);

void Python_finalize(std::mutex* Eigen_model_mutex, std::mutex* training_data_buffer_mutex,
                     std::condition_variable* training_data_buffer_full, bool* start_training, bool* end_training);

int Python_Keras_load_model(std::string model_file_path);

std::vector<double> Python_Keras_predict(std::vector<std::vector<double>> x, int batch_size);  

void training_data_buffer_append(std::vector<std::vector<double>>& training_data_buffer,
                                 std::vector<std::vector<double>>& new_values);

int Python_Keras_training_thread(EigenModel* Eigen_model,
                                 std::mutex* Eigen_model_mutex,
                                 TrainingData* training_data_buffer,
                                 std::mutex* training_data_buffer_mutex,
                                 std::condition_variable* training_data_buffer_full,
                                 bool* start_training, bool* end_training,
                                 const RuntimeParameters& params);

void Python_Keras_set_weights_as_Eigen(EigenModel& eigen_model);

Eigen::MatrixXd eigen_inference_batched(const Eigen::Ref<Eigen::MatrixXd>& input_batch, const EigenModel& model);

std::vector<double> Eigen_predict(const EigenModel& model, std::vector<std::vector<double>> x, int batch_size,
                                  std::mutex* Eigen_model_mutex);

// Otherwise, define the necessary stubs
#else
inline void Python_Keras_setup(std::string functions_file_path){}
inline void set_ai_surrogate_runtime_params(RInsidePOET&, RuntimeParameters&, InitialList&){}
inline void Python_finalize(std::mutex*, std::mutex*, std::condition_variable*, bool*, bool*){}
inline void Python_Keras_load_model(std::string model_file_path){}
inline void training_data_buffer_append(std::vector<std::vector<double>>&, std::vector<std::vector<double>>&){}
inline std::vector<double> Python_Keras_predict(std::vector<std::vector<double>>, int){return {};}
inline int Python_Keras_training_thread(EigenModel*, std::mutex*, 
                                        TrainingData*, std::mutex*,
                                        std::condition_variable*, RuntimeParameters&,
                                        bool*, bool*){return {};}
inline void Python_Keras_set_weights_as_Eigen(EigenModel&){}
inline std::vector<double> Eigen_predict(const EigenModel&, std::vector<std::vector<double>>, int, std::mutex*){return {};}
#endif
} // namespace poet

#endif // AI_FUNCTIONS_HPP