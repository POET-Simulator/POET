/*
** Copyright (C) 2018-2021 Alexander Lindemann, Max Luebke (University of
** Potsdam)
**
** Copyright (C) 2018-2022 Marco De Lucia, Max Luebke (GFZ Potsdam)
**
** POET is free software; you can redistribute it and/or modify it under the
** terms of the GNU General Public License as published by the Free Software
** Foundation; either version 2 of the License, or (at your option) any later
** version.
**
** POET is distributed in the hope that it will be useful, but WITHOUT ANY
** WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
** A PARTICULAR PURPOSE. See the GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License along with
** this program; if not, write to the Free Software Foundation, Inc., 51
** Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include "poet/DHT_Types.hpp"
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <poet/SimParams.hpp>

#include <RInside.h>
#include <Rcpp.h>

#include <iostream>
#include <string>
#include <string_view>
#include <vector>

using namespace poet;
using namespace std;
using namespace Rcpp;

poet::GridParams::s_GridParams(RInside &R) {
  auto tmp_n_cells =
      Rcpp::as<std::vector<uint32_t>>(R.parseEval("mysetup$grid$n_cells"));
  assert(tmp_n_cells.size() < 3);

  this->dim = tmp_n_cells.size();

  std::copy(tmp_n_cells.begin(), tmp_n_cells.end(), this->n_cells.begin());

  auto tmp_s_cells =
      Rcpp::as<std::vector<double>>(R.parseEval("mysetup$grid$s_cells"));

  assert(tmp_s_cells.size() == this->dim);

  std::copy(tmp_s_cells.begin(), tmp_s_cells.end(), this->s_cells.begin());

  this->total_n =
      (dim == 1 ? this->n_cells[0] : this->n_cells[0] * this->n_cells[1]);

  this->type = Rcpp::as<std::string>(R.parseEval("mysetup$grid$type"));
  this->init_df =
      Rcpp::as<Rcpp::DataFrame>(R.parseEval("mysetup$grid$init_cell"));
  this->props =
      Rcpp::as<std::vector<std::string>>(R.parseEval("mysetup$grid$props"));
  this->input_script =
      Rcpp::as<std::string>(R.parseEval("mysetup$grid$input_script"));
  this->database_path =
      Rcpp::as<std::string>(R.parseEval("mysetup$grid$database"));
}

poet::DiffusionParams::s_DiffusionParams(RInside &R) {
  this->initial_t =
      Rcpp::as<Rcpp::DataFrame>(R.parseEval("mysetup$diffusion$init"));
  this->alpha =
      Rcpp::as<Rcpp::NumericVector>(R.parseEval("mysetup$diffusion$alpha"));
  if (Rcpp::as<bool>(
          R.parseEval("'vecinj_inner' %in% names(mysetup$diffusion)"))) {
    this->vecinj_inner =
        Rcpp::as<Rcpp::List>(R.parseEval("mysetup$diffusion$vecinj_inner"));
  }
  this->vecinj =
      Rcpp::as<Rcpp::DataFrame>(R.parseEval("mysetup$diffusion$vecinj"));
  this->vecinj_index =
      Rcpp::as<Rcpp::DataFrame>(R.parseEval("mysetup$diffusion$vecinj_index"));
}

void poet::ChemistryParams::initFromR(RInsidePOET &R) {
  this->database_path =
      Rcpp::as<std::string>(R.parseEval("mysetup$chemistry$database"));
  this->input_script =
      Rcpp::as<std::string>(R.parseEval("mysetup$chemistry$input_script"));

  if (R.checkIfExists("dht_species", "mysetup$chemistry")) {
    this->dht_signifs =
        R.wrapNamedVector<std::uint32_t>("mysetup$chemistry$dht_species");
  }

  if (R.checkIfExists("pht_species", "mysetup$chemistry")) {
    this->pht_signifs =
        R.wrapNamedVector<std::uint32_t>("mysetup$chemistry$pht_species");
  }
}

SimParams::SimParams(int world_rank_, int world_size_) {
  this->simparams.world_rank = world_rank_;
  this->simparams.world_size = world_size_;
}

int SimParams::parseFromCmdl(char *argv[], RInsidePOET &R) {
  // initialize argh object
  argh::parser cmdl(argv);

  // if user asked for help
  if (cmdl[{"help", "h"}]) {
    if (simparams.world_rank == 0) {
      cout << "Todo" << endl
           << "See README.md for further information." << endl;
    }
    return poet::PARSER_HELP;
  }
  // if positional arguments are missing
  else if (!cmdl(2)) {
    if (simparams.world_rank == 0) {
      cerr << "ERROR. Kin needs 2 positional arguments: " << endl
           << "1) the R script defining your simulation and" << endl
           << "2) the directory prefix where to save results and profiling"
           << endl;
    }
    return poet::PARSER_ERROR;
  }

  // collect all parameters which are not known, print them to stderr and return
  // with PARSER_ERROR
  std::list<std::string> optionsError = validateOptions(cmdl);
  if (!optionsError.empty()) {
    if (simparams.world_rank == 0) {
      cerr << "Unrecognized option(s):\n" << endl;
      for (auto option : optionsError) {
        cerr << option << endl;
      }
      cerr << "\nMake sure to use available options. Exiting!" << endl;
    }
    return poet::PARSER_ERROR;
  }

  simparams.print_progressbar = cmdl[{"P", "progress"}];

  /*Parse DHT arguments*/
  chem_params.use_dht = cmdl["dht"];
  // cout << "CPP: DHT is " << ( dht_enabled ? "ON" : "OFF" ) << '\n';

  if (chem_params.use_dht) {
    cmdl("dht-size", DHT_SIZE_PER_PROCESS_MB) >> chem_params.dht_size;
    // cout << "CPP: DHT size per process (Byte) = " << dht_size_per_process <<
    // endl;

    cmdl("dht-snaps", 0) >> chem_params.dht_snaps;

    cmdl("dht-file") >> chem_params.dht_file;
  }
  /*Parse work package size*/
  cmdl("work-package-size", WORK_PACKAGE_SIZE_DEFAULT) >> simparams.wp_size;

  chem_params.use_interp = cmdl["interp"];
  cmdl("interp-size", 100) >> chem_params.pht_size;
  cmdl("interp-min", 5) >> chem_params.interp_min_entries;
  cmdl("interp-bucket-entries", 20) >> chem_params.pht_max_entries;

  /*Parse output options*/
  simparams.store_result = !cmdl["ignore-result"];

  if (simparams.world_rank == 0) {
    cout << "CPP: Complete results storage is "
         << (simparams.store_result ? "ON" : "OFF") << endl;
    cout << "CPP: Work Package Size: " << simparams.wp_size << endl;
    cout << "CPP: DHT is " << (chem_params.use_dht ? "ON" : "OFF") << '\n';

    if (chem_params.use_dht) {
      cout << "CPP: DHT strategy is " << simparams.dht_strategy << endl;
      cout << "CPP: DHT key default digits (ignored if 'signif_vector' is "
              "defined) = "
           << simparams.dht_significant_digits << endl;
      cout << "CPP: DHT logarithm before rounding: "
           << (simparams.dht_log ? "ON" : "OFF") << endl;
      cout << "CPP: DHT size per process (Byte) = " << chem_params.dht_size
           << endl;
      cout << "CPP: DHT save snapshots is " << chem_params.dht_snaps << endl;
      cout << "CPP: DHT load file is " << chem_params.dht_file << endl;
    }
  }

  cmdl(1) >> filesim;
  cmdl(2) >> out_dir;

  chem_params.dht_outdir = out_dir;

  /* distribute information to R runtime */
  // if local_rank == 0 then master else worker
  R["local_rank"] = simparams.world_rank;
  // assign a char* (string) to 'filesim'
  R["filesim"] = wrap(filesim);
  // assign a char* (string) to 'fileout'
  R["fileout"] = wrap(out_dir);
  // pass the boolean "store_result" to the R process
  R["store_result"] = simparams.store_result;
  // worker count
  R["n_procs"] = simparams.world_size - 1;
  // work package size
  R["work_package_size"] = simparams.wp_size;
  // dht enabled?
  R["dht_enabled"] = chem_params.use_dht;
  // log before rounding?
  R["dht_log"] = simparams.dht_log;

  // eval the init string, ignoring any returns
  R.parseEvalQ("source(filesim)");
  R.parseEvalQ("mysetup <- setup");

  this->chem_params.initFromR(R);

  return poet::PARSER_OK;
}

void SimParams::initVectorParams(RInside &R) {}

std::list<std::string> SimParams::validateOptions(argh::parser cmdl) {
  /* store all unknown parameters here */
  std::list<std::string> retList;

  /* loop over all flags and compare to given flaglist*/
  for (auto &flag : cmdl.flags()) {
    if (!(flaglist.find(flag) != flaglist.end()))
      retList.push_back(flag);
  }

  /* and loop also over params and compare to given paramlist */
  for (auto &param : cmdl.params()) {
    if (!(paramlist.find(param.first) != paramlist.end()))
      retList.push_back(param.first);
  }

  return retList;
}
