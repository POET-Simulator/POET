#include "InitialList.hpp"

namespace poet {
void InitialList::initializeFromList(const Rcpp::List &setup) {
  initGrid(setup["grid"]);
  initDiffusion(setup["diffusion"]);
}
} // namespace poet