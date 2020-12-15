#ifndef RRUNTIME_H
#define RRUNTIME_H

#include <RInside.h>
#include <Rcpp.h>
#include <string>

namespace poet {

/**
 * RRuntime is a wrapper class around a RInside (R) runtime and provides several
 * simplified methods to use R commands inside POET.
 *
 * If an instance of RRuntime is created a R runtime will also be spawned.
 */
class RRuntime : public RInside {

public:
  /**
   * Constructor of class RRuntime calling constructor of RInside.
   */
  RRuntime(const int argc, const char *const argv[]) : RInside(argc, argv){};

  void to_C_domain(double *buffer);
  void from_C_domain(double *buffer);

  void setBufferDataFrame(std::string dfname);
  Rcpp::DataFrame getBufferDataFrame();
  size_t getBufferNCol();
  size_t getBufferNRow();

private:
  Rcpp::DataFrame dfbuff;
};
} // namespace poet
#endif // RRUNTIME_H
