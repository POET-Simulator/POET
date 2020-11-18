#include "r_utils.h"

/* This function converts a pure double dataframe into a double array.
    buffer      <- double array, needs to be allocated before
    df          <- reference to input dataframe
*/
void convert_R_Dataframe_2_C_buffer(double* buffer, Rcpp::DataFrame &df)
{
    size_t rowCount = df.nrow();
    size_t colCount = df.ncol();

    for (size_t i = 0; i < rowCount; i++)
    {
        for (size_t j = 0; j < colCount; j++)
        {
            /* Access column vector j and extract value of line i */
            Rcpp::DoubleVector col = df[j];
            buffer[i * colCount + j] = col[i];
        }
    }
}

/* This function converts a double array into a double dataframe.
    buffer      <- input double array
    df          <- reference to output dataframe, needs to be of fitting size, structure will be taken from it
*/
void convert_C_buffer_2_R_Dataframe(double* buffer, Rcpp::DataFrame &df)
{
    size_t rowCount = df.nrow();
    size_t colCount = df.ncol();

    for (size_t i = 0; i < rowCount; i++)
    {
        for (size_t j = 0; j < colCount; j++)
        {
            /* Access column vector j and extract value of line i */
            Rcpp::DoubleVector col = df[j];
            col[i] = buffer[i * colCount + j];
        }
    }
}



