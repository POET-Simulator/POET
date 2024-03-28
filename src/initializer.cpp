#include "Init/InitialList.hpp"
#include "poet.hpp"

#include <Rcpp/vector/instantiation.h>
#include <cstdlib>

#include <RInside.h>
#include <Rcpp.h>

int main(int argc, char **argv) {
  if (argc < 2 || argc > 2) {
    std::cerr << "Usage: " << argv[0] << " <script.R>\n";
    return 1;
  }

  RInside R(argc, argv);

  R.parseEvalQ(init_r_library);

  const std::string script = argv[1];
  Rcpp::Function source_R("source");

  source_R(script);

  Rcpp::List setup = R["setup"];

  //   Rcpp::List grid = R.parseEval("setup$grid");

  //   Rcpp::List results;

  poet::InitialList init(R);

  init.initializeFromList(setup);

  Rcpp::Function save("saveRDS");

  save(init.exportList(), "init.rds");

  //   parseGrid(R, grid, results);
  return EXIT_SUCCESS;
}