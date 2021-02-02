#include "Grid.h"

using namespace poet;
using namespace Rcpp;

void Grid::init() {
  R.parseEval("GRID_TMP <- mysetup$state_C");
  this->ncol = R.parseEval("ncol(GRID_TMP)");
  this->nrow = R.parseEval("nrow(GRID_TMP)");
}

unsigned int Grid::getCols() { return this->ncol; }

unsigned int Grid::getRows() { return this->nrow; }

void Grid::shuffleAndExport(double *buffer) {
  R.parseEval("tmp <- shuffle_field(mysetup$state_T, ordered_ids)");
  R.setBufferDataFrame("tmp");
  R.to_C_domain(buffer);
}

void Grid::importAndUnshuffle(double *buffer) {
  R.setBufferDataFrame("GRID_TMP");
  R.from_C_domain(buffer);
  R["GRID_CHEM_DATA"] = R.getBufferDataFrame();
  R.parseEval("result <- unshuffle_field(GRID_CHEM_DATA, ordered_ids)");
}

void Grid::importWP(double *buffer, unsigned int wp_size) {
  R["GRID_WP_SKELETON"] = getSkeletonDataFrame(wp_size);
  R.setBufferDataFrame("GRID_WP_SKELETON");
  R.from_C_domain(buffer);
  R["work_package_full"] = R.getBufferDataFrame();
}
void Grid::exportWP(double *buffer) {
  R.setBufferDataFrame("result_full");
  R.to_C_domain(buffer);
}

Rcpp::DataFrame Grid::getSkeletonDataFrame(unsigned int rows) {
  R["GRID_ROWS"] = rows;

  Rcpp::DataFrame tmp = R.parseEval("head(GRID_TMP,GRID_ROWS)");
  return tmp;
}
