#ifndef RPOET_H_
#define RPOET_H_

#include "DataStructures.hpp"

#include <RInside.h>
#include <Rcpp.h>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

class RInsidePOET : public RInside {
public:
  static RInsidePOET &getInstance() {
    static RInsidePOET instance;

    return instance;
  }

  RInsidePOET(RInsidePOET const &) = delete;
  void operator=(RInsidePOET const &) = delete;

  inline bool checkIfExists(const std::string &R_name,
                            const std::string &where) {
    return Rcpp::as<bool>(
        this->parseEval("'" + R_name + "' %in% names(" + where + ")"));
  }

  template <class T> inline T getSTL(const std::string &R_name) {
    return Rcpp::as<T>(RInside::operator[](R_name));
  }

  template <typename value_type>
  poet::NamedVector<value_type> wrapNamedVector(const std::string &R_name) {
    std::size_t name_size = this->parseEval("length(names(" + R_name + "))");
    if (name_size == 0) {
      throw std::runtime_error("wrapNamedVector: " + R_name +
                               " is not a named vector!");
    }

    auto names = Rcpp::as<std::vector<std::string>>(
        this->parseEval("names(" + R_name + ")"));
    auto values = Rcpp::as<Rcpp::NumericVector>(this->parseEval(R_name));

    poet::NamedVector<value_type> ret_map;

    for (std::size_t i = 0; i < names.size(); i++) {
      ret_map.insert(
          std::make_pair(names[i], static_cast<value_type>(values[i])));
    }

    return ret_map;
  }

private:
  RInsidePOET() : RInside(){};
};

#endif // RPOET_H_
