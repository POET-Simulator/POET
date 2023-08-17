#include "AdvectionModule.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include <Rcpp.h>

namespace poet {

void simulate(double dt);

void AdvectionModule::initializeParams(RInsidePOET &R) {
  const std::uint32_t n = this->n_cells[0] =
      R.parseEval("mysetup$grid$n_cells[1]");
  const std::uint32_t m = this->n_cells[1] =
      R.parseEval("mysetup$grid$n_cells[2]");

  const std::uint32_t field_size{n * m};

  const auto init_vec =
      Rcpp::as<Rcpp::NumericVector>(R.parseEval("mysetup$advection$init"));

  const auto prop_names = Rcpp::as<std::vector<std::string>>(init_vec.names());
  std::vector<std::vector<double>> init_field;
  init_field.reserve(prop_names.size());

  for (const auto &val : init_vec) {
    init_field.push_back(std::vector<double>(field_size, val));
  }

  const Rcpp::DataFrame vecinj = R.parseEval("mysetup$advection$vecinj");
  const Rcpp::DataFrame inner_vecinj =
      R.parseEval("mysetup$advection$vecinj_inner");

  for (std::size_t i = 0; i < inner_vecinj.size(); i++) {
    const Rcpp::NumericVector tuple = inner_vecinj[i];
    const std::uint32_t cell_index = (tuple[2] - 1) + ((tuple[1] - 1) * n);
    for (std::size_t prop_i = 0; prop_i < prop_names.size(); prop_i++) {
      const Rcpp::NumericVector curr_prop_vec_inj = vecinj[prop_names[prop_i]];
      init_field[prop_i][cell_index] = curr_prop_vec_inj[tuple[0] - 1];
    }
    this->cells_const.insert(cell_index);
  }

  this->t_field = Field(field_size, init_field, prop_names);
}

} // namespace poet
