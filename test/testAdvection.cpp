#include "testDataStructures.hpp"
#include <poet/AdvectionModule.hpp>
#include <doctest/doctest.h>

using namespace poet;

TEST_CASE("Advection Module") {
    auto &R = RInsidePOET::getInstance();

    R["input_script"] = SampleInputScript;
    R.parseEvalQ("source(input_script)");
    R.parseEval("mysetup <- setup");

    AdvectionModule adv(R);
}
