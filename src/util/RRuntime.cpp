#include "RRuntime.h"

#include <RInside.h>
#include <Rcpp.h>

#include <string>

using namespace poet;

/**
 * Convert a R dataframe into a C continious memory area.
 *
 * @param varname Name of the R internal variable name.
 */
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

/**
 * Convert continious C memory area into R dataframe and puts it into R runtime.
 *
 * @param buffer Pointer to memory area which should be converted into R
 * dataframe.
 * @param skeleton Defines the raw data frame structure and muste be defined
 * inside the R runtime beforehand.
 * @param varname Name of the R internal variable name.
 */
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
