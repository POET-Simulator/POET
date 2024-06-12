#include <Rcpp.h>
#include <cstddef>
#include <cstdint>
#include <doctest/doctest.h>
#include <string>
#include <vector>

#include <Base/RInsidePOET.hpp>
#include <DataStructures/Field.hpp>

#include "testDataStructures.hpp"

using namespace poet;

#define CHECK_AND_FAIL_LOOP(val1, val2)                                        \
  if (val1 != val2) {                                                          \
    FAIL_CHECK("Value differs: ", val1, " != ", val2);                         \
  }

TEST_CASE("Field") {
  constexpr uint32_t VEC_SIZE = 5;
  constexpr uint32_t VEC_COUNT = 3;
  constexpr double INIT_VAL = 1;

  std::vector<std::string> names = {"C", "Ca", "Na"};
  std::vector<FieldColumn> init_values(names.size(),
                                       FieldColumn(VEC_SIZE, INIT_VAL));

  RInsidePOET &R = RInsidePOET::getInstance();

  R["sourcefile"] = RInside_source_file;
  R.parseEval("source(sourcefile)");

  SUBCASE("Initialize field with 2D vector") {
    Field dut(VEC_SIZE, init_values, names);

    auto props = dut.GetProps();

    CHECK_EQ(names.size(), props.size());

    const auto res = dut["C"];

    CHECK_EQ(res.size(), VEC_SIZE);

    for (const auto &elem : res) {
      CHECK_AND_FAIL_LOOP(elem, INIT_VAL);
    }
  }

  SUBCASE("Initialize field with 2D vector") {
    std::vector<double> init_values(VEC_SIZE * VEC_COUNT, 1);
    Field dut(VEC_SIZE, init_values, names);

    auto props = dut.GetProps();

    CHECK_EQ(names.size(), props.size());

    const auto res = dut["C"];

    CHECK_EQ(res.size(), VEC_SIZE);

    for (const auto &elem : res) {
      CHECK_AND_FAIL_LOOP(elem, INIT_VAL);
    }
  }

  Field dut(VEC_SIZE, init_values, names);

  SUBCASE("Return vector") {
    std::vector<double> output = dut.AsVector();

    CHECK(output.size() == VEC_SIZE * VEC_COUNT);
  }

  SUBCASE("SEXP return") {
    R["test"] = dut.asSEXP();

    Field field_constructed = R.parseEval("test");

    Rcpp::DataFrame R_df = R.parseEval("test");
    Rcpp::CharacterVector R_names = R_df.names();

    for (std::size_t i = 0; i < VEC_COUNT; i++) {
      const auto r_values = Rcpp::as<std::vector<double>>(R_df[i]);
      CHECK_EQ(r_values, dut[names[i]]);
      CHECK_EQ(R_names[i], names[i]);
      CHECK_EQ(field_constructed[names[i]], dut[names[i]]);
    }
  }

  SUBCASE("Apply R function (set Na to zero)") {
    poet::DEFunc to_call("simple_field");
    Field field_proc = to_call(dut.asSEXP());

    CHECK_EQ(field_proc["Na"], FieldColumn(dut.GetRequestedVecSize(), 0));
  }

  SUBCASE("Apply R function (add two fields)") {
    poet::DEFunc to_call("extended_field");
    Field field_proc = to_call(dut.asSEXP(), dut.asSEXP());

    CHECK_EQ(field_proc["Na"],
             FieldColumn(dut.GetRequestedVecSize(), INIT_VAL + INIT_VAL));
  }

  constexpr double NEW_VAL = 2.;
  std::vector<double> new_val_vec(VEC_SIZE, NEW_VAL);

  dut["C"] = new_val_vec;

  SUBCASE("Check manipulation of column") {

    auto test_it = new_val_vec.begin();

    for (const auto &val : dut["C"]) {
      CHECK_EQ(val, *test_it);
      test_it++;
    }
  }

  SUBCASE("Check correctly manipulated values of 1D vector") {
    auto out_res = dut.AsVector();
    auto out_it = out_res.begin();

    for (uint32_t i = 0; i < VEC_SIZE; i++, out_it++) {
      CHECK_AND_FAIL_LOOP(NEW_VAL, *out_it);
    }

    for (; out_it != out_res.end(); out_it++) {
      CHECK_AND_FAIL_LOOP(INIT_VAL, *out_it);
    }
  }

  std::vector<double> new_field(VEC_SIZE * VEC_COUNT);

  for (uint32_t i = 0; i < VEC_COUNT; i++) {
    for (uint32_t j = 0; j < VEC_SIZE; j++) {
      new_field[j + (i * VEC_SIZE)] = (double)(i + 1) / (double)(j + 1);
    }
  }

  dut = new_field;

  SUBCASE("SetFromVector return correct field vector") {
    auto ret_vec = dut.AsVector();

    auto ret_it = ret_vec.begin();
    auto new_it = new_field.begin();

    for (; ret_it != ret_vec.end(); ret_it++, new_it++) {
      CHECK_AND_FAIL_LOOP(*ret_it, *new_it);
    }
  }

  SUBCASE("Get single column with new values") {
    auto new_it = new_field.begin() + 2 * VEC_SIZE;
    auto ret_vec = dut["Na"];
    auto ret_it = ret_vec.begin();

    CHECK_EQ(ret_vec.size(), VEC_SIZE);

    for (; ret_it != ret_vec.end(); ret_it++, new_it++) {
      CHECK_AND_FAIL_LOOP(*ret_it, *new_it);
    }
  }

  std::vector<double> mg_vec = {1, 2, 3, 4, 5};
  dut["Mg"] = mg_vec;

  SUBCASE("Operator creates new element and places it at the end") {

    auto new_props = dut.GetProps();

    CHECK_EQ(new_props.size(), 4);
    CHECK_EQ(new_props[3], "Mg");

    auto ret_vec = dut.As2DVector();
    auto mg_vec_it = mg_vec.begin();

    for (const auto &val : ret_vec[3]) {
      CHECK_AND_FAIL_LOOP(val, *(mg_vec_it++));
    }
  }

  // reset field
  names = dut.GetProps();

  dut = std::vector<FieldColumn>(names.size(), FieldColumn(VEC_SIZE, INIT_VAL));

  constexpr double SOME_OTHER_VAL = -0.5;
  std::vector<std::string> some_other_props = {"Na", "Cl"};
  std::vector<double> some_other_values(VEC_SIZE * some_other_props.size(),
                                        SOME_OTHER_VAL);

  Field some_other_field(VEC_SIZE, some_other_values, some_other_props);

  SUBCASE("Update existing field from another field") {
    dut.update(some_other_field);

    auto ret_vec = dut.As2DVector();
    auto ret_it = ret_vec.begin();

    for (const auto &prop : names) {
      const auto &curr_vec = *(ret_it++);

      for (const auto &val : curr_vec) {
        CHECK_AND_FAIL_LOOP((prop == "Na" ? SOME_OTHER_VAL : INIT_VAL), val);
      }
    }
  }
}
