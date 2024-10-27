/**
 * @file AI_functions.hpp
 * @author Hans Straile (straile@uni-potsdam.de)
 * @brief API for the AI/Machine Learning based chemistry surrogate model with functions to initialize a neural network and use it for training and inference via Keras for Python . 
 * @version 0.1
 * @date 01 Nov 2024
 *
 * This file implements functions to train and predict with a neural network.
 * All functions are based on a user supplied Keras model. Pyhton.h is used for model
 * training with Keras and can be used for inference. The default inference funtion 
 * is implemented with Eigen matrices in C++. All functions use 2 different models
 * to process data separately according to a K-means cluster assignement. This file
 * alo contains the functions for the K-means algorithm.
 *  
 */

#ifndef AI_FUNCTIONS_H
#define AI_FUNCTIONS_H

#include <string>
#include <vector>
#include "poet.hpp"

// PhreeqC definition of pi clashes with Eigen macros
// so we have to temporarily undef it 
#pragma push_macro("pi")
#undef pi
#include <Eigen/Dense> 
#pragma pop_macro("pi")

namespace poet {

struct EigenModel {
  // The first model will be used for all values if clustering is disabled
  // or for the reactive part of the field if clustering is enabled
  std::vector<Eigen::MatrixXd> weight_matrices;
  std::vector<Eigen::VectorXd> biases;

  // The other model will be used for the non-reactive cluster
  // (if clustering is enabled)
  std::vector<Eigen::MatrixXd> weight_matrices_no_reaction;
  std::vector<Eigen::VectorXd> biases_no_reaction;
};

struct TrainingData {
  std::vector<std::vector<double>> x;
  std::vector<std::vector<double>> y;
  std::vector<int> cluster_labels;
  int n_training_runs = 0;
};

// Ony declare the actual functions if flag is set 
#ifdef USE_AI_SURROGATE

int Python_Keras_setup(std::string functions_file_path, std::string cuda_src_dir);

void Python_finalize(std::mutex* Eigen_model_mutex, std::mutex* training_data_buffer_mutex,
                     std::condition_variable* training_data_buffer_full, bool* start_training, 
                     bool* end_training);

int Python_Keras_load_model(std::string model, std::string model_reactive,
                            bool use_clustering);

std::vector<int> K_Means(std::vector<std::vector<double>>& field, int k, int maxIterations = 100);

std::vector<double> Python_Keras_predict(std::vector<std::vector<double>>& x, int batch_size,
                                         std::vector<int>& cluster_labels);  

void training_data_buffer_append(std::vector<std::vector<double>>& training_data_buffer,
                                 std::vector<std::vector<double>>& new_values);

void cluster_labels_append(std::vector<int>& labels, std::vector<int>& new_labels,
                           std::vector<int> validity);

int Python_Keras_training_thread(EigenModel* Eigen_model, EigenModel* Eigen_model_reactive,
                                 std::mutex* Eigen_model_mutex,
                                 TrainingData* training_data_buffer,
                                 std::mutex* training_data_buffer_mutex,
                                 std::condition_variable* training_data_buffer_full,
                                 bool* start_training, bool* end_training,
                                 const RuntimeParameters& params);

void update_weights(EigenModel* model, const std::vector<std::vector<std::vector<double>>>& weights);

std::vector<std::vector<std::vector<double>>> Python_Keras_get_weights(std::string model_name);

std::vector<double> Eigen_predict_clustered(const EigenModel& model, const EigenModel& model_reactive,
                                            std::vector<std::vector<double>>& x, 
                                            int batch_size, std::mutex* Eigen_model_mutex,
                                            std::vector<int>& cluster_labels);
std::vector<double> Eigen_predict(const EigenModel& model, std::vector<std::vector<double>> x, int batch_size,
                                  std::mutex* Eigen_model_mutex);


// Otherwise, define the necessary stubs
#else
inline void Python_Keras_setup(std::string, std::string){}
inline void Python_finalize(std::mutex*, std::mutex*, std::condition_variable*, bool*, bool*){}
inline void Python_Keras_load_model(std::string, std::string, bool){}
inline std::vector<int> K_Means(std::vector<std::vector<double>>&, int, int) {return {};}
inline std::vector<double> Python_Keras_predict(std::vector<std::vector<double>>&, int,
                                                std::vector<int>&){return {};}
inline void training_data_buffer_append(std::vector<std::vector<double>>&,
                                        std::vector<std::vector<double>>&){}
inline void cluster_labels_append(std::vector<int>&, std::vector<int>&, std::vector<int>){}
inline int Python_Keras_training_thread(EigenModel*, , EigenModel*, std::mutex*, 
                                        TrainingData*, std::mutex*, std::condition_variable*,
                                        bool*, bool*, const RuntimeParameters&){return {};}

inline void update_weights(EigenModel*, const std::vector<std::vector<std::vector<double>>>&){}
inline std::vector<std::vector<std::vector<double>>> Python_Keras_get_weights(std::string){return {};}
inline std::vector<double> Eigen_predict_clustered(const EigenModel&, const EigenModel&,
                                                   std::vector<std::vector<double>>&, int,
                                                   std::mutex*, std::vector<int>&){return {};}
inline std::vector<double> Eigen_predict(const EigenModel&, std::vector<std::vector<double>>, int,
                                         std::mutex*){return {};}
#endif
} // namespace poet

#endif // AI_FUNCTIONS_HPP