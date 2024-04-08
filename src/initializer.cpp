#include "Init/InitialList.hpp"
#include "poet.hpp"

#include <Rcpp/vector/instantiation.h>
#include <cstdlib>

#include <RInside.h>
#include <Rcpp.h>
#include <iostream>

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

  poet::InitialList init(R);

  init.initializeFromList(setup);

  // replace file extension by .rds
  const std::string rds_out_filename =
      script.substr(0, script.find_last_of('.')) + ".rds";

  Rcpp::Function save("saveRDS");
  save(init.exportList(), Rcpp::wrap(rds_out_filename));

  std::cout << "Saved result to " << rds_out_filename << std::endl;

  //   parseGrid(R, grid, results);
  return EXIT_SUCCESS;
}