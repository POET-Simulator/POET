#include "InitialList.hpp"

namespace poet {
void InitialList::initializeFromList(const Rcpp::List &setup) {
  initGrid(setup[grid_key]);
  initDiffusion(setup[diffusion_key]);
}
} // namespace poet