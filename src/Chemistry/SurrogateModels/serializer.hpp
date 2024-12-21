#ifndef SERIALIZER_H
#define SERIALIZER_H

#include "AI_functions.hpp"
#include <cstddef>

namespace poet{

/**
 * @brief Serialize the weights and biases of the model into a memory location
 * to send them via RDMA
 *
 * @param model: Struct of EigenModel containing the weights and biases of the
 * model
 * @param memory: Pointer to the memory location where the serialized data will
 * be stored
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
 * @return std::vector<char> 
 */
int serializeTrainingData(const TrainingData& data, void *memory);

/**
 * @brief Deserialize the training data from a memory location
 * 
 * @param data Pointer to the memory location where the serialized data is stored
 * @return TrainingData struct containing the training data
 */
TrainingData deserializeTrainingData(void* data);

/**
 * @brief Calculates the size of stored elements in the EigenModel and TrainData
 * structs
 *
 * @param struct_pointer: pointer to the struct
 * @param type: determines the struct type given: E for EigenModel and T TrainData
 * @return size_t: size of stored elements
 */
size_t calculateStructSize(void* struct_pointer, char type);
}
#endif