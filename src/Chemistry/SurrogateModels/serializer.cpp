#include "serializer.hpp"
#include "AI_functions.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <vector>

using namespace std;
namespace poet{
size_t calculateStructSize(void *struct_pointer, char type){

    size_t struct_size = 0;
    if (type == 'E') {
      struct_size += sizeof(size_t); // number of matrices
      struct_size +=
          static_cast<EigenModel *>(struct_pointer)->weight_matrices.size() *
          2 * sizeof(size_t);        // dimensions of matrices
      struct_size += sizeof(size_t); // number of vectors
      struct_size += static_cast<EigenModel *>(struct_pointer)->biases.size() *
                     sizeof(size_t); // length of vectors
      
      for (const Eigen::MatrixXd &matrix :
           static_cast<EigenModel *>(struct_pointer)->weight_matrices) {
        // fprintf(stderr, "matrix size: rows:%td, cols: %td\n", matrix.rows(), matrix.cols());
        struct_size += matrix.size() * sizeof(double);
        // fprintf(stderr, "matrix size %td\n", matrix.size());

      }
      for (const Eigen::VectorXd &bias :
           static_cast<EigenModel *>(struct_pointer)->biases) {
        struct_size += bias.size() * sizeof(double);
        // fprintf(stderr, "matrix size %td\n", bias.size());

      }
    } else if (type == 'T') {
        struct_size += sizeof(size_t); // number of vectors
        struct_size += sizeof(size_t); // length of vector
        for (const std::vector<double> &vector:
            *static_cast<std::vector<std::vector<double>>*>(struct_pointer)){
                struct_size += vector.size() * sizeof(double);
            }
    }

return struct_size;

}


int serializeModelWeights(const EigenModel *model, char *memory){
   
    size_t num_matrices = model->weight_matrices.size();
    size_t size_counter = 0;
    std::memcpy(memory, &num_matrices, sizeof(size_t));
    memory += sizeof(size_t);
    size_counter += sizeof(size_t);
    for (const Eigen::MatrixXd &matrix : model->weight_matrices) {
      size_t rows = matrix.rows(), cols = matrix.cols();
      fprintf(stdout, "rows: %zu, cols: %zu\n", rows, cols);
      std::memcpy(memory, &rows, sizeof(size_t));
      memory += sizeof(size_t);
      size_counter += sizeof(size_t);
      std::memcpy(memory, &cols, sizeof(size_t));
      memory += sizeof(size_t);
      size_counter += sizeof(size_t);
      std::memcpy(memory, matrix.data(), rows * cols * sizeof(double));
      memory += rows * cols * sizeof(double);
      size_counter += rows * cols * sizeof(double);
    }

    // Serialisierung der Bias-Vektoren
    size_t num_biases = model->biases.size();
    std::memcpy(memory, &num_biases, sizeof(size_t));
    memory += sizeof(size_t);
    size_counter += sizeof(size_t);
    for (const Eigen::VectorXd &bias : model->biases) {
      size_t size = bias.size();
      std::memcpy(memory, &size, sizeof(size_t));
      memory += sizeof(size_t);
      size_counter += sizeof(size_t);
      std::memcpy(memory, bias.data(), size * sizeof(double));
      memory += size * sizeof(double);
      size_counter += size * sizeof(double);
    }
    fprintf(stdout, "serializer size: %zu\n", size_counter);
    return 0;
}

EigenModel deserializeModelWeights(char *memory, size_t buffer_size) {
    EigenModel deserializedModel;
    size_t size_counter = 0;

    // Anzahl Matrizen
    size_t num_matrices;
    if (buffer_size < sizeof(size_t)) {
        fprintf(stderr, "Buffer too small.\n");
        return deserializedModel;
    }
    std::memcpy(&num_matrices, memory, sizeof(size_t));
    memory += sizeof(size_t);
    size_counter += sizeof(size_t);
    deserializedModel.weight_matrices.resize(num_matrices);

    // matrix deserialization
    for (Eigen::MatrixXd &matrix : deserializedModel.weight_matrices) {
        size_t rows, cols;

        // buffer check
        if (size_counter + 2 * sizeof(size_t) > buffer_size) {
            fprintf(stderr, "Buffer too small for matrix dimensions.\n");
            return deserializedModel;
        }

        std::memcpy(&rows, memory, sizeof(size_t));
        memory += sizeof(size_t);
        size_counter += sizeof(size_t);

        std::memcpy(&cols, memory, sizeof(size_t));
        memory += sizeof(size_t);
        size_counter += sizeof(size_t);

        if (size_counter + rows * cols * sizeof(double) > buffer_size) {
            fprintf(stderr, "Buffer too small for matrix data.\n");
            return deserializedModel;
        }

        // interpret memory as Eigen::MatrixXd (more efficient than memcpy?)
        Eigen::MatrixXd temp = Eigen::Map<Eigen::MatrixXd>(
            reinterpret_cast<double*>(memory), rows, cols);
        matrix = temp;  // copy data to matrix
        memory += rows * cols * sizeof(double);
        size_counter += rows * cols * sizeof(double);
    }

    // number of bias vectors
    size_t num_biases;
    if (size_counter + sizeof(size_t) > buffer_size) {
        fprintf(stderr, "Buffer too small for biases.\n");
        return deserializedModel;
    }
    std::memcpy(&num_biases, memory, sizeof(size_t));
    memory += sizeof(size_t);
    size_counter += sizeof(size_t);
    deserializedModel.biases.resize(num_biases);

    // deserialization of bias vectors
    for (Eigen::VectorXd &bias : deserializedModel.biases) {
        size_t size;
        if (size_counter + sizeof(size_t) > buffer_size) {
            fprintf(stderr, "Buffer too small for bias size.\n");
            return deserializedModel;
        }

        std::memcpy(&size, memory, sizeof(size_t));
        memory += sizeof(size_t);
        size_counter += sizeof(size_t);

        if (size_counter + size * sizeof(double) > buffer_size) {
            fprintf(stderr, "Buffer too small for bias data.\n");
            return deserializedModel;
        }

        // same procedure as for the matrices
        // TODO: delete temp variable
        Eigen::VectorXd temp = Eigen::Map<Eigen::VectorXd>(
            reinterpret_cast<double*>(memory), size);
        bias = temp;  // Kopieren der Daten
        memory += size * sizeof(double);
        size_counter += size * sizeof(double);
    }

    return deserializedModel;
}


int serializeTrainingData(std::vector<std::vector<double>> data, char *memory){

    size_t num_vectors = data.size();

    std::memcpy(memory, &num_vectors, sizeof(size_t));
    memory += sizeof(size_t);

    for (const std::vector<double> &vector : data) {
        size_t size = vector.size();
        std::memcpy(memory, &size, sizeof(size_t));
        memory += sizeof(size_t);
        std::memcpy(memory, vector.data(), size * sizeof(double));
        memory += size * sizeof(double);
    }

    return 0;
}

std::vector<std::vector<double>> deserializeTrainingData(char* data){

    std::vector<std::vector<double>> deserialized_data;
    size_t num_vectors;
    std::memcpy(&num_vectors, data, sizeof(size_t));
    data += sizeof(size_t);

    for (size_t i = 0; i < num_vectors; i++) {
        size_t size;
        std::memcpy(&size, data, sizeof(size_t));
        data += sizeof(size_t);

        std::vector<double> vector(size);
        std::memcpy(vector.data(), data, size * sizeof(double));
        data += size * sizeof(double);

        deserialized_data.push_back(vector);
    }

    return deserialized_data;

}

}