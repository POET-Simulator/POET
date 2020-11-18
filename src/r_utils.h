#pragma once
#include <RInside.h> 

/*Functions*/
void convert_R_Dataframe_2_C_buffer(double* buffer, Rcpp::DataFrame &df);
void convert_C_buffer_2_R_Dataframe(double* buffer, Rcpp::DataFrame &df);