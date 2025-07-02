#pragma once 

#include <cstdint>
#include <DataStructures/Field.hpp>

struct Checkpoint_s {
    poet::Field &field;
    uint32_t iteration;
};