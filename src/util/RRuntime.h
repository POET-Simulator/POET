#ifndef RRUNTIME_H
#define RRUNTIME_H

#include <RInside.h>
#include <Rcpp.h>

#include <string>

namespace poet {

/**
 * @brief Provides an interface to a R runtime.
 *
 * RRuntime is a wrapper class around a RInside (R) runtime and provides several
 * simplified methods to use R commands inside POET.
 *
 */
class RRuntime : public RInside {
 public:
  /**
   * @brief Construct a new RRuntime object
   *
   * Since this is an inherited class of RInside the constructor of the super
   * class is just called.
   *
   * @param argc Argument counter of the program
   * @param argv Argument values of the program
   */
  RRuntime(const int argc, const char *const argv[]) : RInside(argc, argv){};

  /**
   * Convert a R dataframe into a C continious memory area.
   *
   * @param buffer Name of the R internal variable name.
   */

  /**
   * @brief Convert a R dataframe into a C continious memory area.
   *
   * A buffer data frame must be set beforehand with setBufferDataFrame. Then
   * each value will be set into the continious memory area. This is done row
   * wise.
   *
   * @todo: Might be more performant if all columns would be loaded at once and
   * not for each column seperatly
   *
   * @param buffer Pointer to pre-allocated memory
   */
  void to_C_domain(double *buffer);

  /**
   * @brief Convert continious C memory area into R dataframe and puts it into R
   * runtime.
   *
   * A buffer data frame must be set beforehand with setBufferDataFrame. Then
   * each value will be set into buffered data frame of this object. This is
   * done row wise.
   *
   * @todo: Might be more performant if all columns would be loaded at once and
   * not for each column seperatly
   *
   * @param buffer Pointer to memory area which should be converted into R
   * dataframe.
   */
  void from_C_domain(double *buffer);

  /**
   * @brief Set the Buffer Data Frame object
   *
   * Set the buffered data frame (will be mostly the grid) of this object.
   *
   * @param dfname Name of the data frame inside R runtime
   */
  void setBufferDataFrame(std::string dfname);

  /**
   * @brief Get the Buffer Data Frame object
   *
   * Returning the current buffered data frame as a Rcpp data frame.
   *
   * @return Rcpp::DataFrame Current buffered data frame
   */
  Rcpp::DataFrame getBufferDataFrame();

  /**
   * @brief Get the Buffer N Col object
   *
   * Get the numbers of columns of the buffered data frame.
   *
   * @return size_t Count of columns of buffered data frame
   */
  size_t getBufferNCol();

  /**
   * @brief Get the Buffer N Row object
   *
   * Get the numbers of rows of the buffered data frame.
   *
   * @return size_t Count of rows of buffered data frame
   */
  size_t getBufferNRow();

 private:
  /**
   * @brief Buffered data frame
   *
   * This is used to convert a R data frame into continious memory used by C/C++
   * runtime and vice versa. Must be set with setBufferDataFrame and can be
   * manipulated with from_C_domain.
   * 
   * @todo: Find a cleaner solution. Maybe abstraction via another class.
   *
   */
  Rcpp::DataFrame dfbuff;
};
}  // namespace poet
#endif  // RRUNTIME_H
