#include <doctest/doctest.h>

#include <Chemistry/SurrogateModels/Rounding.hpp>

TEST_CASE("Rounding") {
  constexpr int signif = 3;

  poet::DHT_Rounder rounder;

  SUBCASE("+exp - +sig") {
    double input = 1.2345;
    auto rounded = rounder.round(input, signif, false);
    CHECK_EQ(rounded.sc_notation.significant, 123);
    CHECK_EQ(rounded.sc_notation.exp, 0);
  }

  SUBCASE("+exp - -sig") {
    double input = -1.2345;
    auto rounded = rounder.round(input, signif, false);
    CHECK_EQ(rounded.sc_notation.significant, -123);
    CHECK_EQ(rounded.sc_notation.exp, 0);
  }

  SUBCASE("-exp - +sig") {
    double input = 1.23456789E-5;
    auto rounded = rounder.round(input, signif, false);
    CHECK_EQ(rounded.sc_notation.significant, 123);
    CHECK_EQ(rounded.sc_notation.exp, -5);
  }

  SUBCASE("-exp - -sig") {
    double input = -1.23456789E-5;
    auto rounded = rounder.round(input, signif, false);
    CHECK_EQ(rounded.sc_notation.significant, -123);
    CHECK_EQ(rounded.sc_notation.exp, -5);
  }

  SUBCASE("-exp - +sig (exceeding aqueous minimum)") {
    double input = 1.23456789E-14;
    auto rounded = rounder.round(input, signif, false);
    CHECK_EQ(rounded.sc_notation.significant, 0);
    CHECK_EQ(rounded.sc_notation.exp, 0);
  }

  SUBCASE("-exp - -sig (exceeding aqueous minimum)") {
    double input = -1.23456789E-14;
    auto rounded = rounder.round(input, signif, false);
    CHECK_EQ(rounded.sc_notation.significant, 0);
    CHECK_EQ(rounded.sc_notation.exp, 0);
  }

  SUBCASE("-exp - +sig (diff exceeds aqueous minimum)") {
    double input = 1.23456789E-13;
    auto rounded = rounder.round(input, signif, false);
    CHECK_EQ(rounded.sc_notation.significant, 1);
    CHECK_EQ(rounded.sc_notation.exp, -13);
  }

  SUBCASE("-exp - -sig (diff exceeds aqueous minimum)") {
    double input = -1.23456789E-13;
    auto rounded = rounder.round(input, signif, false);
    CHECK_EQ(rounded.sc_notation.significant, -1);
    CHECK_EQ(rounded.sc_notation.exp, -13);
  }

  SUBCASE("-exp - +sig (diff exceeds aqueous minimum - but H or O)") {
    double input = 1.23456789E-13;
    auto rounded = rounder.round(input, signif, true);
    CHECK_EQ(rounded.sc_notation.significant, 123);
    CHECK_EQ(rounded.sc_notation.exp, -13);
  }

  SUBCASE("-exp - -sig (diff exceeds aqueous minimum - but H or O)") {
    double input = -1.23456789E-13;
    auto rounded = rounder.round(input, signif, true);
    CHECK_EQ(rounded.sc_notation.significant, -123);
    CHECK_EQ(rounded.sc_notation.exp, -13);
  }
}
