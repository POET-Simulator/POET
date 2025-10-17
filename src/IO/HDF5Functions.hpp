#pragma once

#include <string>
#include "Datatypes.hpp"


int write_checkpoint(const std::string &dir_path, const std::string &file_name, struct Checkpoint_s &&checkpoint);

int read_checkpoint(const std::string &dir_path, const std::string &file_name, struct Checkpoint_s &checkpoint);
