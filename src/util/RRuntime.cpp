/*
** Copyright (C) 2018-2021 Alexander Lindemann, Max Luebke (University of
** Potsdam)
**
** Copyright (C) 2018-2021 Marco De Lucia (GFZ Potsdam)
**
** POET is free software; you can redistribute it and/or modify it under the
** terms of the GNU General Public License as published by the Free Software
** Foundation; either version 2 of the License, or (at your option) any later
** version.
**
** POET is distributed in the hope that it will be useful, but WITHOUT ANY
** WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
** A PARTICULAR PURPOSE. See the GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License along with
** this program; if not, write to the Free Software Foundation, Inc., 51
** Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

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
