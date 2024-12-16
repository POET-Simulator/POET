#include "Init/InitialList.hpp"
#include "poet.hpp"

#include <CLI/CLI.hpp>

#include <cstdlib>

#include <RInside.h>
#include <Rcpp.h>
#include <iostream>
#include <string>

int main(int argc, char **argv) {
  // initialize RIinside
  RInside R(argc, argv);

  R.parseEvalQ(init_r_library);
  R.parseEvalQ(kin_r_library);

  // parse command line arguments
  CLI::App app{"POET Initializer - Poet R script to POET qs/rds converter"};

  std::string input_script;
  app.add_option("PoetScript.R", input_script, "POET R script to convert")
      ->required();

  std::string output_file;
  app.add_option("-n, --name", output_file,
                 "Name of the output file without file extension");

  bool setwd;
  app.add_flag("-s, --setwd", setwd,
               "Set working directory to the directory of the directory of the "
               "input script")
      ->default_val(false);

  bool asRDS{false};
  app.add_flag("-r, --rds", asRDS, "Save output as .rds")->default_val(false);

  bool asQS{false};
  app.add_flag("-q, --qs", asQS, "Save output as .qs")->default_val(false);

  bool includeH0O0{false};
  app.add_flag("--include-h0-o0", includeH0O0,
               "Include H(0) and O(0) in the output")
      ->default_val(false);

  CLI11_PARSE(app, argc, argv);

  // source the input script
  std::string normalized_path_script;
  std::string in_file_name;

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

  // if output file is not specified, use the input file name
  if (output_file.empty()) {
    std::string curr_path =
        Rcpp::as<std::string>(Rcpp::Function("normalizePath")(Rcpp::wrap(".")));

    output_file = curr_path + "/" +
                  in_file_name.substr(0, in_file_name.find_last_of('.'));
  }

  // append the correct file extension
  if (asRDS) {
    output_file += ".rds";
  } else if (asQS) {
    output_file += ".qs";
  } else {
    output_file += ".qs2";
  }

  // set working directory to the directory of the input script
  if (setwd) {
    const std::string dir_path = Rcpp::as<std::string>(
        Rcpp::Function("dirname")(normalized_path_script));

    Rcpp::Function("setwd")(Rcpp::wrap(dir_path));
  }

  Rcpp::List setup = R["setup"];

  poet::InitialList init(R);

  init.initializeFromList(setup);

  // use the generic handler defined in kin_r_library.R
  Rcpp::Function SaveRObj_R("SaveRObj");
  SaveRObj_R(init.exportList(), Rcpp::wrap(output_file));

  std::cout << "Saved result to " << output_file << std::endl;

  //   parseGrid(R, grid, results);
  return EXIT_SUCCESS;
}
