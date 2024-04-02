#ifndef RPOET_H_
#define RPOET_H_

#include <RInside.h>
#include <Rcpp.h>
#include <Rinternals.h>
#include <cstddef>
#include <exception>
#include <optional>
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

private:
  RInsidePOET() : RInside(){};
};

template <typename T> class RHookFunction {
public:
  RHookFunction() {}
  RHookFunction(RInside &R, const std::string &f_name) {
    try {
      this->func = Rcpp::Function(Rcpp::as<SEXP>(R.parseEval(f_name.c_str())));
    } catch (const std::exception &e) {
    }
  }

  RHookFunction(SEXP f) {
    try {
      this->func = Rcpp::Function(f);
    } catch (const std::exception &e) {
    }
  }

  template <typename... Args> T operator()(Args... args) const {
    if (func.has_value()) {
      return (Rcpp::as<T>(this->func.value()(args...)));
    } else {
      throw std::exception();
    }
  }

  bool isValid() const { return this->func.has_value(); }

  SEXP asSEXP() const { return Rcpp::as<SEXP>(this->func.value()); }

private:
  std::optional<Rcpp::Function> func;
};

#endif // RPOET_H_
