#pragma once

#include <Rcpp.h>

namespace poet {
template <typename T> class NamedVector : public Rcpp::NumericVector {
public:
  NamedVector() : Rcpp::NumericVector(){};

  NamedVector(const std::vector<std::string> &in_names,
              const std::vector<T> &in_values)
      : Rcpp::NumericVector(Rcpp::wrap(in_values)) {
    this->names() = Rcpp::CharacterVector(Rcpp::wrap(in_names));
  }

  NamedVector(const SEXP &s) : Rcpp::NumericVector(s){};

  bool empty() const { return (this->size() == 0); }

  std::vector<std::string> getNames() const {
    return Rcpp::as<std::vector<std::string>>(this->names());
  }
  std::vector<T> getValues() const { return Rcpp::as<std::vector<T>>(*this); }
};
} // namespace poet