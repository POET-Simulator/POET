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

// PhreeqC definition of pi clashes with Eigen macros so we have to temporarily undef it 
#pragma push_macro("pi")
#undef pi
#include <Eigen/Dense> 
#pragma pop_macro("pi")

namespace poet {

struct EigenModel {
    std::vector<Eigen::MatrixXd> weight_matrices;
    std::vector<Eigen::VectorXd> biases;
};

struct TrainingData {
  std::vector<std::vector<double>> x;
  std::vector<std::vector<double>> y;
  int n_training_runs = 0;
};

// Ony declare the actual functions if flag is set 
#ifdef USE_AI_SURROGATE

int Python_Keras_setup(std::string functions_file_path, std::string cuda_src_dir);

void Python_finalize(std::mutex* Eigen_model_mutex, std::mutex* training_data_buffer_mutex,
                     std::condition_variable* training_data_buffer_full, bool* start_training, bool* end_training);

int Python_Keras_load_model(std::string model_reaction, std::string model_no_reaction);

std::vector<int> kMeans(std::vector<std::vector<double>>& field, int k, int maxIterations = 100);

std::vector<double> Python_Keras_predict(std::vector<std::vector<double>>& x, int batch_size);  

void training_data_buffer_append(std::vector<std::vector<double>>& training_data_buffer,
                                 std::vector<std::vector<double>>& new_values);

int Python_Keras_training_thread(EigenModel* Eigen_model,
                                 std::mutex* Eigen_model_mutex,
                                 TrainingData* training_data_buffer,
                                 std::mutex* training_data_buffer_mutex,
                                 std::condition_variable* training_data_buffer_full,
                                 bool* start_training, bool* end_training,
                                 const RuntimeParameters& params);

void update_weights(EigenModel* model, const std::vector<std::vector<std::vector<double>>>& weights);

std::vector<std::vector<std::vector<double>>> Python_Keras_get_weights();

std::vector<double> Eigen_predict(const EigenModel& model, std::vector<std::vector<double>>& x, int batch_size,
                                  std::mutex* Eigen_model_mutex);

// Otherwise, define the necessary stubs
#else
inline void Python_Keras_setup(std::string, std::string){}
inline void Python_finalize(std::mutex*, std::mutex*, std::condition_variable*, bool*, bool*){}
inline void Python_Keras_load_model(std::string, std::string){}
inline std::vector<int> kMeans(std::vector<std::vector<double>>&, int, int) {return {};}
inline std::vector<double> Python_Keras_predict(std::vector<std::vector<double>>&, int){return {};}
inline void training_data_buffer_append(std::vector<std::vector<double>>&, std::vector<std::vector<double>>&){}
inline int Python_Keras_training_thread(EigenModel*, std::mutex*, 
                                        TrainingData*, std::mutex*,
                                        std::condition_variable*,
                                        bool*, bool*, const RuntimeParameters&){return {};}

inline void update_weights(EigenModel*, const std::vector<std::vector<std::vector<double>>>&){}
inline std::vector<std::vector<std::vector<double>>> Python_Keras_get_weights(){return {};}
inline std::vector<double> Eigen_predict(const EigenModel&, std::vector<std::vector<double>>&, int, std::mutex*){return {};}
#endif
} // namespace poet

#endif // AI_FUNCTIONS_HPP