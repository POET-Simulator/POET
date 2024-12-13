
#ifndef LOOKUPKEY_H_
#define LOOKUPKEY_H_

#include "HashFunctions.hpp"

#include <cstdint>
#include <cstring>
#include <vector>

namespace poet {

constexpr std::int8_t SC_NOTATION_EXPONENT_MASK = -128;
constexpr std::int64_t SC_NOTATION_SIGNIFICANT_MASK = 0xFFFFFFFFFFFF;

struct Lookup_SC_notation {
  std::int8_t exp : 8;
  std::int64_t significant : 56;

  constexpr static Lookup_SC_notation nan() {
    return {SC_NOTATION_EXPONENT_MASK, SC_NOTATION_SIGNIFICANT_MASK};
  }

  constexpr bool isnan() {
    return !!(exp == SC_NOTATION_EXPONENT_MASK &&
              significant == SC_NOTATION_SIGNIFICANT_MASK);
  }
};

union Lookup_Keyelement {
  double fp_element;
  Lookup_SC_notation sc_notation;

  bool operator==(const Lookup_Keyelement &other) const {
    return std::memcmp(this, &other, sizeof(Lookup_Keyelement)) == 0 ? true
                                                                     : false;
  }
};

class LookupKey : public std::vector<Lookup_Keyelement> {
public:
  using std::vector<Lookup_Keyelement>::vector;

  std::vector<double> to_double() const;
  static Lookup_SC_notation round_from_double(const double in,
                                              std::uint32_t signif);
  static double to_double(const Lookup_SC_notation in);
};

struct LookupKeyHasher {
  std::uint64_t operator()(const LookupKey &k) const {
    const uint32_t key_size = k.size() * sizeof(Lookup_Keyelement);

    return poet::Murmur2_64A(key_size, k.data());
  }
};

} // namespace poet

#endif // LOOKUPKEY_H_
