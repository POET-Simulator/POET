#include "Init/InitialList.hpp"
#include "poet.hpp"

#include "Base/argh.hpp"

#include <cstdlib>

#include <RInside.h>
#include <Rcpp.h>
#include <iostream>
#include <string>

int main(int argc, char **argv) {

  argh::parser cmdl({"-o", "--output"});

  cmdl.parse(argc, argv);

  if (cmdl[{"-h", "--help"}] || cmdl.pos_args().size() != 2) {
    std::cout << "Usage: " << argv[0]
              << " [-o, --output output_file]"
                 " [-s, --setwd] "
                 " <script.R>"
              << std::endl;
    return EXIT_SUCCESS;
  }

  RInside R(argc, argv);

  R.parseEvalQ(init_r_library);

  std::string input_script = cmdl.pos_args()[1];
  std::string normalized_path_script;
  std::string in_file_name;

  std::string curr_path =
      Rcpp::as<std::string>(Rcpp::Function("normalizePath")(Rcpp::wrap(".")));

  try {
    normalized_path_script =
        Rcpp::as<std::string>(Rcpp::Function("normalizePath")(input_script));

    in_file_name = Rcpp::as<std::string>(
        Rcpp::Function("basename")(normalized_path_script));

    std::cout << "Using script " << input_script << std::endl;
    Rcpp::Function("source")(input_script);
  } catch (Rcpp::exception &e) {
    std::cerr << "Error while sourcing " << input_script << e.what()
              << std::endl;
    return EXIT_FAILURE;
  }

  std::string output_file;

  cmdl({"-o", "--output"},
       curr_path + "/" +
           in_file_name.substr(0, in_file_name.find_last_of('.')) + ".rds") >>
      output_file;

  if (cmdl[{"-s", "--setwd"}]) {
    const std::string dir_path = Rcpp::as<std::string>(
        Rcpp::Function("dirname")(normalized_path_script));

    Rcpp::Function("setwd")(Rcpp::wrap(dir_path));
  }

  Rcpp::List setup = R["setup"];

  poet::InitialList init(R);

  init.initializeFromList(setup);

  // replace file extension by .rds
  Rcpp::Function save("saveRDS");

  save(init.exportList(), Rcpp::wrap(output_file));

  std::cout << "Saved result to " << output_file << std::endl;

  //   parseGrid(R, grid, results);
  return EXIT_SUCCESS;
}