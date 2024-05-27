#pragma once

#include <cstdint>
#include <vector>

namespace poet {

enum DHT_PROP_TYPES { DHT_TYPE_DEFAULT, DHT_TYPE_CHARGE, DHT_TYPE_TOTAL };
enum CHEMISTRY_OUT_SOURCE { CHEM_PQC, CHEM_DHT, CHEM_INTERP, CHEM_AISURR };

struct WorkPackage {
  std::size_t size;
  std::vector<std::vector<double>> input;
  std::vector<std::vector<double>> output;
  std::vector<std::uint8_t> mapping;

  WorkPackage(std::size_t _size) : size(_size) {
    input.resize(size);
    output.resize(size);
    mapping.resize(size, CHEM_PQC);
  }
};
} // namespace poet