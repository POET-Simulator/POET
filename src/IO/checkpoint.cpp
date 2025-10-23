#include "IO/Datatypes.hpp"
#include <cstdint>
#include <highfive/H5Easy.hpp>

int write_checkpoint(const std::string &file_path, struct Checkpoint_s &&checkpoint){

    // TODO: errorhandling
    H5Easy::File file(file_path, H5Easy::File::Overwrite);


    H5Easy::dump(file, "/MetaParam/Iterations", checkpoint.iteration);
    H5Easy::dump(file, "/Grid/Names", checkpoint.field.GetProps());
    H5Easy::dump(file, "/Grid/Chemistry", checkpoint.field.As2DVector());

    return 0;
}

int read_checkpoint(const std::string &file_path, struct Checkpoint_s &checkpoint){

     H5Easy::File file(file_path, H5Easy::File::ReadOnly);

     checkpoint.iteration = H5Easy::load<uint32_t>(file, "/MetaParam/Iterations");

     checkpoint.field = H5Easy::load<std::vector<std::vector<double>>>(file, "/Grid/Chemistry");

    return 0;
}

