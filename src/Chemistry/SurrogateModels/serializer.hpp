#ifndef SERIALIZER_H
#define SERIALIZER_H

#include "AI_functions.hpp"
#include <cstddef>

namespace poet{

/**
 * @brief Calculates the size of stored elements in the EigenModel and TrainData
 * structs
 *
 * @param struct_pointer: pointer to the struct
 * @param type: determines the struct type given: E for EigenModel and T
 * training data vector structures
 * @return size_t: size of stored elements
 */
size_t calculateStructSize(void* struct_pointer, char type);

/**
 * @brief Serialize the weights and biases of the model into a memory location
 * to send them via RDMA
 *
 * @param model: Struct of EigenModel containing the weights and biases of the
 * model
 * @param memory: Pointer to the memory location where the serialized data will
 * be stored
 * The serialized data looks like this:
 * |# matrices|# cols of matrix 1|# rows of matrix 1|matrix 1 data|# cols of matrix 2|...
 * |# bias vectors|length of bias vector 1|bias vector 1 data|length of bias vector 2|...
 * @return int: 0 if the serialization was successful, -1 otherwise
 */
int serializeModelWeights(const EigenModel *model, char *memory);

/**
 * @brief Deserialize the weights and biases of the model from a memory location
 * 
 * @param data Pointer to the memory location where the serialized data is stored
 * @return EigenModel struct containing the weights and biases of the model
 */
EigenModel deserializeModelWeights(char* memory, size_t buffer_size);

/**
 * @brief Serialize the training data into a memory location to send it via RDMA
 * 
 * @param data Struct of TrainingData containing the training data
 * @param memory
 * The serialized data looks like this:
 * |# of vectors|length of vector 1.1|vector 1.1 data|length of vector 1.2|...
 * |# of vectors|length of vector 2.1|vector 2.1 data|length of vector 2.2|...
 * |length of vector |vector 3 data|
 * n_training_runs
 * @return std::vector<char> 
 */
int serializeTrainingData(std::vector<std::vector<double>> *data, char *memory);

/**
 * @brief Deserialize the training data from a memory location
 *
 * @param data Pointer to the memory location where the serialized data is
 * stored
 * @return std::vector<std::vector<double>> containing n vectors for each
 * species with m training elements
 */
std::vector<std::vector<double>> deserializeTrainingData(char* data);

/**
 * @brief Serialize the weights and biases of the model into a memory location
 * to send them via RDMA
 *
 * @param model: 3d vector containing the weights and biases of the
 * model
 * @param memory: Pointer to the memory location where the serialized data will
 * be stored
 * The serialized data looks as follows:
 * |# layers|# rows of matrix 1|# cols of matrix 1|matrix 1 data|# rows of matrix 2|...
 * @return int: 0 if the serialization was successful, -1 otherwise
 */
int serializeCPPWeights(std::vector<std::vector<std::vector<double>>> &cpp_weights, char* memory);

std::vector<std::vector<std::vector<double>>> deserializeCPPWeights(char* data);
}
#endif