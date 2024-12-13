#ifndef ROUNDING_H_
#define ROUNDING_H_

#include "LookupKey.hpp"
#include <cmath>
#include <cstdint>

namespace poet {

constexpr std::int8_t AQUEOUS_EXP = -13;

template <typename Input, typename Output, typename ConvertTo = double>
class IRounding {
public:
  virtual Output round(const Input &, std::uint32_t signif) = 0;
  virtual ConvertTo convert(const Output &) = 0;
};

class DHT_Rounder {
public:
  Lookup_Keyelement round(const double &value, std::uint32_t signif,
                          bool is_ho) {

    if (std::isnan(value)) {
      return {.sc_notation = Lookup_SC_notation::nan()};
    }

    std::int8_t exp =
        static_cast<std::int8_t>(std::floor(std::log10(std::fabs(value))));

    if (!is_ho) {
      if (exp < AQUEOUS_EXP) {
        return {.sc_notation = {0, 0}};
      }

      std::int8_t diff = exp - signif + 1;
      if (diff < AQUEOUS_EXP) {
        signif -= AQUEOUS_EXP - diff;
      }
    }

    Lookup_Keyelement elem;
    elem.sc_notation.exp = exp;
    elem.sc_notation.significant =
        static_cast<std::int64_t>(value * std::pow(10, signif - exp - 1));

    return elem;
  }

  double convert(const Lookup_Keyelement &key_elem) {
    std::int32_t normalized_exp = static_cast<std::int32_t>(
        -std::log10(std::fabs(key_elem.sc_notation.significant)));

    // add stored exponent to normalized exponent
    normalized_exp += key_elem.sc_notation.exp;

    // return significant times 10 to the power of exponent
    return key_elem.sc_notation.significant * std::pow(10., normalized_exp);
  }
};

class PHT_Rounder : public IRounding<Lookup_Keyelement, Lookup_Keyelement> {
public:
  Lookup_Keyelement round(const Lookup_Keyelement &value,
                          std::uint32_t signif) {
    Lookup_Keyelement new_val = value;

    std::uint32_t diff_signif =
        static_cast<std::uint32_t>(
            std::ceil(std::log10(std::abs(value.sc_notation.significant)))) -
        signif;

    new_val.sc_notation.significant = static_cast<int64_t>(
        value.sc_notation.significant / std::pow(10., diff_signif));

    if (new_val.sc_notation.significant == 0) {
      new_val.sc_notation.exp = 0;
    }

    return new_val;
  }

  double convert(const Lookup_Keyelement &key_elem) {
    std::int32_t normalized_exp = static_cast<std::int32_t>(
        -std::log10(std::fabs(key_elem.sc_notation.significant)));

    // add stored exponent to normalized exponent
    normalized_exp += key_elem.sc_notation.exp;

    // return significant times 10 to the power of exponent
    return key_elem.sc_notation.significant * std::pow(10., normalized_exp);
  }
};

} // namespace poet

#endif // ROUNDING_H_
