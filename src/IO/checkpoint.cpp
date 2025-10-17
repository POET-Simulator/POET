#include "IO/Datatypes.hpp"
#include <cstdint>
#include <highfive/H5Easy.hpp>
#include <filesystem>

namespace fs = std::filesystem;

int write_checkpoint(const std::string &dir_path, const std::string &file_name, struct Checkpoint_s &&checkpoint){
    
    if (!fs::exists(dir_path)) {
        std::cerr << "Directory does not exist: " << dir_path << std::endl;
        return -1;
    }
    fs::path file_path = fs::path(dir_path) / file_name;
    // TODO: errorhandling
    H5Easy::File file(file_path, H5Easy::File::Overwrite);

    H5Easy::dump(file, "/MetaParam/Iterations", checkpoint.iteration);
    H5Easy::dump(file, "/Grid/Names", checkpoint.field.GetProps());
    H5Easy::dump(file, "/Grid/Chemistry", checkpoint.field.As2DVector());

    return 0;
}

int read_checkpoint(const std::string &dir_path, const std::string &file_name, struct Checkpoint_s &checkpoint){
    
    fs::path file_path = fs::path(dir_path) / file_name;

    if (!fs::exists(file_path)) {
        std::cerr << "File does not exist: " << file_path << std::endl;
        return -1;
    }
     
    H5Easy::File file(file_path, H5Easy::File::ReadOnly);

     checkpoint.iteration = H5Easy::load<uint32_t>(file, "/MetaParam/Iterations");

     checkpoint.field = H5Easy::load<std::vector<std::vector<double>>>(file, "/Grid/Chemistry");

    return 0;
}

