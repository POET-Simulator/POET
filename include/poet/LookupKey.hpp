//  Time-stamp: "Last modified 2023-08-11 10:12:52 mluebke"

#ifndef LOOKUPKEY_H_
#define LOOKUPKEY_H_

#include "poet/HashFunctions.hpp"
#include <cstdint>
#include <cstring>
#include <vector>

namespace poet {

struct Lookup_SC_notation {
  std::int8_t exp : 8;
  std::int64_t significant : 56;
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
