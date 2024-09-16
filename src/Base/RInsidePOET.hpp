#pragma once

#include <RInside.h>
#include <Rcpp.h>
#include <Rinternals.h>
#include <exception>
#include <memory>
#include <string>

namespace poet {
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

/**
 * @brief Deferred evaluation function
 *
 * The class is intended to call R functions within an existing RInside
 * instance. The problem with "original" Rcpp::Function is that they require:
 * 1. RInside instance already present, restricting the declaration of
 *     Rcpp::Functions in global scope
 * 2. Require the function to be present. Otherwise, they will throw an
 *     exception.
 * This class solves both problems by deferring the evaluation of the function
 *  until the constructor is called and evaluating whether the function is
 *  present or not, wihout throwing an exception.
 *
 * @tparam T Return type of the function
 */
class DEFunc {
public:
  DEFunc() {}
  DEFunc(const std::string &f_name) {
    try {
      this->func = std::make_shared<Rcpp::Function>(f_name);
    } catch (const std::exception &e) {
    }
  }

  DEFunc(SEXP f) {
    try {
      this->func = std::make_shared<Rcpp::Function>(f);
    } catch (const std::exception &e) {
    }
  }

  template <typename... Args> SEXP operator()(Args... args) const {
    if (func) {
      return (*this->func)(args...);
    } else {
      throw std::exception();
    }
  }

  DEFunc &operator=(const DEFunc &rhs) {
    this->func = rhs.func;
    return *this;
  }

  DEFunc(const DEFunc &rhs) { this->func = rhs.func; }

  bool isValid() const { return static_cast<bool>(func); }

  SEXP asSEXP() const {
    if (!func) {
      return R_NilValue;
    }
    return Rcpp::as<SEXP>(*this->func.get());
  }

private:
  std::shared_ptr<Rcpp::Function> func;
};

} // namespace poet