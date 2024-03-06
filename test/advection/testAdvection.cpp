#include <Transport/AdvectionModule.hpp>

#include <cstddef>
#include <string>

#include "InputFiles.hpp"

using namespace poet;

constexpr std::size_t MAX_ITER = 10;
constexpr double DT = 200;

int main(int argc, char **argv) {
  auto &R = RInsidePOET::getInstance();

  R["input_script"] = SampleInputScript;
  R.parseEvalQ("source(input_script)");
  R.parseEval("mysetup <- setup");

  AdvectionModule adv(R);
  R["field"] = adv.getField().asSEXP();
  R.parseEval("saveRDS(field, 'adv_0.rds')");

  for (std::size_t i = 0; i < MAX_ITER; i++) {
    adv.simulate(DT);
    const std::string save_q =
        "saveRDS(field, 'adv_" + std::to_string(i + 1) + ".rds')";
    R["field"] = adv.getField().asSEXP();
    R.parseEval(save_q);
  }
}
