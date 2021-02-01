#include "RRuntime.h"

#include <RInside.h>
#include <Rcpp.h>

#include <string>

using namespace poet;

void RRuntime::to_C_domain(double *buffer) {
  size_t rowCount = dfbuff.nrow();
  size_t colCount = dfbuff.ncol();

  for (size_t i = 0; i < rowCount; i++) {
    for (size_t j = 0; j < colCount; j++) {
      /* Access column vector j and extract value of line i */
      Rcpp::DoubleVector col = dfbuff[j];
      buffer[i * colCount + j] = col[i];
    }
  }
}

void RRuntime::from_C_domain(double *buffer) {
  size_t rowCount = dfbuff.nrow();
  size_t colCount = dfbuff.ncol();

  for (size_t i = 0; i < rowCount; i++) {
    for (size_t j = 0; j < colCount; j++) {
      /* Access column vector j and extract value of line i */
      Rcpp::DoubleVector col = dfbuff[j];
      col[i] = buffer[i * colCount + j];
    }
  }
}

void RRuntime::setBufferDataFrame(std::string dfname) {
  this->dfbuff = parseEval(dfname);
}

Rcpp::DataFrame RRuntime::getBufferDataFrame() { return this->dfbuff; }

size_t RRuntime::getBufferNCol() { return (this->dfbuff).ncol(); }

size_t RRuntime::getBufferNRow() { return (this->dfbuff).nrow(); }
