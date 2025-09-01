#pragma once

#include <string>
#include "Datatypes.hpp"

int write_checkpoint(const std::string &file_path, struct Checkpoint_s &&checkpoint);

int read_checkpoint(const std::string &file_path, struct Checkpoint_s &checkpoint);

