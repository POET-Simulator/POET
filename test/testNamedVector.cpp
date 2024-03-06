#include <Rcpp.h>
#include <cstddef>
#include <doctest/doctest.h>
#include <utility>
#include <vector>

#include <Base/RInsidePOET.hpp>
#include <DataStructures/DataStructures.hpp>

#include <InputFiles.hpp>

TEST_CASE("NamedVector") {
  RInsidePOET &R = RInsidePOET::getInstance();

  R["sourcefile"] = RInside_source_file;
  R.parseEval("source(sourcefile)");

  const std::vector<std::string> names{{"H", "O", "Ca"}};
  const std::vector<double> values{{0.1, 0.2, 0.3}};

  poet::NamedVector<double> nv(names, values);

  SUBCASE("SEXP conversion") {
    R["test"] = nv;

    const Rcpp::NumericVector R_value = R.parseEval("test");
    const Rcpp::CharacterVector R_names = R_value.names();

    const poet::NamedVector<double> nv_constructed = R["test"];

    for (std::size_t i = 0; i < R_value.size(); i++) {
      CHECK_EQ(R_value[i], values[i]);
      CHECK_EQ(R_names[i], names[i]);
      CHECK_EQ(nv_constructed[i], values[i]);
      CHECK_EQ(nv_constructed.getNames()[i], names[i]);
    }
  }

  SUBCASE("Apply R function (set to zero)") {
    RHookFunction<poet::NamedVector<double>> to_call(R, "simple_named_vec");
    nv = to_call(nv);

    CHECK_EQ(nv[2], 0);
  }

  SUBCASE("Apply R function (second NamedVector)") {
    RHookFunction<poet::NamedVector<double>> to_call(R, "extended_named_vec");

    const std::vector<std::string> names{{"C", "H", "Mg"}};
    const std::vector<double> values{{0, 1, 2}};

    poet::NamedVector<double> second_nv(names, values);

    nv = to_call(nv, second_nv);

    CHECK_EQ(nv[0], 1.1);
  }

  SUBCASE("Apply R function (check if zero)") {
    RHookFunction<bool> to_call(R, "bool_named_vec");

    CHECK_FALSE(to_call(nv));
  }
}
