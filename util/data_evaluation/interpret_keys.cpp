/* This file is intended to be used as input file for 'Rcpp::sourceCPP()'.
 *
 * The provided function will translate our key data structure back into human
 * readable double values also interpretable by R or other languages.
 */

#include <Rcpp.h>

#include <cmath>
#include <cstdint>
#include <vector>

using DHT_Keyelement = struct keyelem {
  std::int8_t exp : 8;
  std::int64_t significant : 56;
};

using namespace Rcpp;

// [[Rcpp::export]]
std::vector<double> rcpp_key_convert(std::vector<double> input) {

  std::vector<double> output;
  output.reserve(input.size());

  for (const double &value : input) {
    DHT_Keyelement currKeyelement = *((DHT_Keyelement *)&value);

    double normalize =
        ((std::int32_t)-std::log10(std::fabs(currKeyelement.significant))) +
        currKeyelement.exp;

    output.push_back(currKeyelement.significant * std::pow(10., normalize));
  }

  return output;
}
