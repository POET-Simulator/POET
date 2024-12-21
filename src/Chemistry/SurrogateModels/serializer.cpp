#include "serializer.hpp"
#include "AI_functions.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <cstddef>
#include <cstdio>

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

// EigenModel deserializeModelWeights(char *memory, size_t buffer_size){

//     EigenModel deserializedModel;

    
//     size_t num_matrices; 
//     size_t size_counter = 0;
//     std::memcpy(&num_matrices, memory, sizeof(size_t));
//     fprintf(stdout, "number of matrices: %zu\n", num_matrices);

//     memory += sizeof(size_t);
//     size_counter += sizeof(size_t);
//     deserializedModel.weight_matrices.resize(num_matrices);
//     for (Eigen::MatrixXd &matrix : deserializedModel.weight_matrices) {
//         size_t rows, cols;
//         std::memcpy(&rows, memory, sizeof(size_t));
//         memory += sizeof(size_t);
//         size_counter += sizeof(size_t);
//         std::memcpy(&cols, memory, sizeof(size_t));
//         fprintf(stdout, "rows: %zu, cols: %zu\n", rows, cols);
//         memory += sizeof(size_t);
//         size_counter += sizeof(size_t);
//         fprintf(stdout, "rows before: %td, cols before: %td\n", matrix.rows(), matrix.cols());
//         matrix.resize(rows, cols);
//         std::memcpy(matrix.data(), memory, rows * cols * sizeof(double));
        
//         memory += rows * cols * sizeof(double);
//         size_counter += rows * cols * sizeof(double);
//     }
//     fprintf(stdout, "deserialized size of matrices: %zu\n", size_counter);
//     size_t num_biases;
//     std::memcpy(&num_biases, memory, sizeof(size_t));
//     memory += sizeof(size_t);
//     size_counter += sizeof(size_t);

//     fprintf(stdout, "number of biases: %zu\n", num_biases);
//     deserializedModel.biases.resize(num_biases);
//     for (Eigen::VectorXd &bias : deserializedModel.biases) {
//         size_t size;
//         std::memcpy(&size, memory, sizeof(size_t));
//         fprintf(stdout, "bias length: %zu\n", size);
//         memory += sizeof(size_t);
//         size_counter += sizeof(size_t);
//         bias.resize(size);
//         std::memcpy(bias.data(), memory, size * sizeof(double));
//         memory += size * sizeof(double);
//         size_counter += size * sizeof(double);
//     }
//     fprintf(stdout, "deserialized size: %zu\n", size_counter);
//     if(size_counter > buffer_size){
//       fprintf(stderr, "buffer corrupted!\n");
//     }
//     return deserializedModel;
// }
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

    // Matrizen deserialisieren
    for (Eigen::MatrixXd &matrix : deserializedModel.weight_matrices) {
        size_t rows, cols;

        // Buffer-Check
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

        // Kopiere Daten in neue Matrix
        Eigen::MatrixXd temp = Eigen::Map<Eigen::MatrixXd>(
            reinterpret_cast<double*>(memory), rows, cols);
        matrix = temp;  // Kopieren der Daten
        memory += rows * cols * sizeof(double);
        size_counter += rows * cols * sizeof(double);
    }

    // Anzahl Biases
    size_t num_biases;
    if (size_counter + sizeof(size_t) > buffer_size) {
        fprintf(stderr, "Buffer too small for biases.\n");
        return deserializedModel;
    }
    std::memcpy(&num_biases, memory, sizeof(size_t));
    memory += sizeof(size_t);
    size_counter += sizeof(size_t);
    deserializedModel.biases.resize(num_biases);

    // Biases deserialisieren
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

        // Kopiere Daten in neuen Vektor
        Eigen::VectorXd temp = Eigen::Map<Eigen::VectorXd>(
            reinterpret_cast<double*>(memory), size);
        bias = temp;  // Kopieren der Daten
        memory += size * sizeof(double);
        size_counter += size * sizeof(double);
    }

    return deserializedModel;
}


}