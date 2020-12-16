#ifndef GRID_H
#define GRID_H

#include "../util/RRuntime.h"
#include <Rcpp.h>

namespace poet {
class Grid {

public:
  Grid(RRuntime &R) : R(R){};
  
  void init();
  
  unsigned int getCols();
  unsigned int getRows();
  
  void shuffleAndExport(double *buffer);
  void importAndUnshuffle(double *buffer);

  void importWP(double *buffer, unsigned int p_size);
  void exportWP(double *buffer);

private:
  RRuntime R;
  unsigned int ncol;
  unsigned int nrow;

  Rcpp::DataFrame getSkeletonDataFrame(unsigned int rows);
};
} // namespace poet
#endif // GRID_H
