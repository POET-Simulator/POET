#pragma once

#include "../DataStructures/Field.hpp"

#include <RInside.h>
#include <Rcpp.h>
#include <Rcpp/vector/instantiation.h>
#include <cstdint>

namespace poet {

class InitialList {
public:
  InitialList(RInside &R) : R(R){};

  void initializeFromList(const Rcpp::List &setup);

  void importList(const std::string &file_name);
  Rcpp::List exportList(const std::string &file_name);

private:
  RInside &R;

  // Grid members
  void initGrid(const Rcpp::List &grid_input);
  std::uint8_t dim;

  std::uint32_t n_cols;
  std::uint32_t n_rows;

  double s_cols;
  double s_rows;

  std::vector<std::uint32_t> constant_cells;

  Rcpp::List initial_grid;

  // No export
  Rcpp::NumericMatrix phreeqc_mat;

  // Initialized by grid, modified by chemistry
  std::map<int, std::string> pqc_raw_dumps;
};
} // namespace poet