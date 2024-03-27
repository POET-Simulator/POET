/*
** Copyright (C) 2018-2021 Alexander Lindemann, Max Luebke (University of
** Potsdam)
**
** Copyright (C) 2018-2022 Marco De Lucia, Max Luebke (GFZ Potsdam)
**
** Copyright (C) 2018-2022 Marco De Lucia (GFZ Potsdam), Max Luebke (University
** of Potsdam)
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

#include "SimParams.hpp"

#include "Base/Macros.hpp"

#include <cstdint>

#include <cstdint>
#include <string>
#include <vector>

#include "Base/RInsidePOET.hpp"
#include "argh.hpp" // Argument handler https://github.com/adishavit/argh

using namespace poet;

RuntimeParameters::RuntimeParameters(RInsidePOET &R, char *argv[],
                                     int world_rank_)
    : world_rank(world_rank_) {

  // initialize argh object
  argh::parser cmdl(argv);

  // if user asked for help
  if (cmdl[{"help", "h"}]) {
    if (this->world_rank == 0) {
      MSG("Todo");
      MSG("See README.md for further information.");
    }

    return poet::PARSER_HELP;
  }
  // if positional arguments are missing
  else if (!cmdl(2)) {
    if (simparams.world_rank == 0) {
      ERRMSG("Kin needs 2 positional arguments: ");
      ERRMSG("1) the R script defining your simulation and");
      ERRMSG("2) the directory prefix where to save results and profiling");
    }
    return poet::PARSER_ERROR;
  }

  // collect all parameters which are not known, print them to stderr and return
  // with PARSER_ERROR
  std::list<std::string> optionsError = validateOptions(cmdl);
  if (!optionsError.empty()) {
    if (simparams.world_rank == 0) {
      ERRMSG("Unrecognized option(s):\n");
      for (auto option : optionsError) {
        ERRMSG(std::string(option));
      }
      ERRMSG("Make sure to use available options. Exiting!");
    }
    return poet::PARSER_ERROR;
  }

  simparams.print_progressbar = cmdl[{"P", "progress"}];

  // simparams.print_progressbar = cmdl[{"P", "progress"}];

  /* Parse DHT arguments */
  chem_params.use_dht = cmdl["dht"];
  chem_params.use_interp = cmdl["interp"];
  // cout << "CPP: DHT is " << ( dht_enabled ? "ON" : "OFF" ) << '\n';

  cmdl("dht-size", DHT_SIZE_PER_PROCESS_MB) >> chem_params.dht_size;
  // cout << "CPP: DHT size per process (Byte) = " << dht_size_per_process <<
  // endl;

  cmdl("dht-snaps", 0) >> chem_params.dht_snaps;

  cmdl("dht-file") >> chem_params.dht_file;
  /*Parse work package size*/
  cmdl("work-package-size", WORK_PACKAGE_SIZE_DEFAULT) >> simparams.wp_size;

  cmdl("interp-size", 100) >> chem_params.pht_size;
  cmdl("interp-min", 5) >> chem_params.interp_min_entries;
  cmdl("interp-bucket-entries", 20) >> chem_params.pht_max_entries;

  /*Parse output options*/
  simparams.store_result = !cmdl["ignore-result"];

  /*Parse work package size*/
  cmdl("work-package-size", WORK_PACKAGE_SIZE_DEFAULT) >> simparams.wp_size;

  chem_params.use_interp = cmdl["interp"];
  cmdl("interp-size", 100) >> chem_params.pht_size;
  cmdl("interp-min", 5) >> chem_params.interp_min_entries;
  cmdl("interp-bucket-entries", 20) >> chem_params.pht_max_entries;

  /*Parse output options*/
  simparams.store_result = !cmdl["ignore-result"];

  if (simparams.world_rank == 0) {
    MSG("Complete results storage is " + BOOL_PRINT(simparams.store_result));
    MSG("Work Package Size: " + std::to_string(simparams.wp_size));
    MSG("DHT is " + BOOL_PRINT(chem_params.use_dht));

    if (chem_params.use_dht) {
      MSG("DHT strategy is " + std::to_string(simparams.dht_strategy));
      // MDL: these should be outdated (?)
      // MSG("DHT key default digits (ignored if 'signif_vector' is "
      // 	"defined) = "
      // 	 << simparams.dht_significant_digits);
      // MSG("DHT logarithm before rounding: "
      // 	 << (simparams.dht_log ? "ON" : "OFF"));
      MSG("DHT size per process (Megabyte) = " +
          std::to_string(chem_params.dht_size));
      MSG("DHT save snapshots is " + BOOL_PRINT(chem_params.dht_snaps));
      MSG("DHT load file is " + chem_params.dht_file);
    }

    if (chem_params.use_interp) {
      MSG("PHT interpolation enabled: " + BOOL_PRINT(chem_params.use_interp));
      MSG("PHT interp-size = " + std::to_string(chem_params.pht_size));
      MSG("PHT interp-min  = " +
          std::to_string(chem_params.interp_min_entries));
      MSG("PHT interp-bucket-entries = " +
          std::to_string(chem_params.pht_max_entries));
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
};

// poet::GridParams::s_GridParams(RInside &R) {
//   auto tmp_n_cells =
//       Rcpp::as<std::vector<uint32_t>>(R.parseEval("mysetup$grid$n_cells"));
//   assert(tmp_n_cells.size() < 3);

//   this->dim = tmp_n_cells.size();

//   std::copy(tmp_n_cells.begin(), tmp_n_cells.end(), this->n_cells.begin());

//   auto tmp_s_cells =
//       Rcpp::as<std::vector<double>>(R.parseEval("mysetup$grid$s_cells"));

//   assert(tmp_s_cells.size() == this->dim);

//   std::copy(tmp_s_cells.begin(), tmp_s_cells.end(), this->s_cells.begin());

//   this->total_n =
//       (dim == 1 ? this->n_cells[0] : this->n_cells[0] * this->n_cells[1]);

//   this->type = Rcpp::as<std::string>(R.parseEval("mysetup$grid$type"));
// }

// poet::DiffusionParams::s_DiffusionParams(RInside &R) {
//   this->initial_t =
//       Rcpp::as<Rcpp::DataFrame>(R.parseEval("mysetup$diffusion$init"));
//   this->alpha =
//       Rcpp::as<Rcpp::NumericVector>(R.parseEval("mysetup$diffusion$alpha"));
//   if (Rcpp::as<bool>(
//           R.parseEval("'vecinj_inner' %in% names(mysetup$diffusion)"))) {
//     this->vecinj_inner =
//         Rcpp::as<Rcpp::List>(R.parseEval("mysetup$diffusion$vecinj_inner"));
//   }
//   this->vecinj =
//       Rcpp::as<Rcpp::DataFrame>(R.parseEval("mysetup$diffusion$vecinj"));
//   this->vecinj_index =
//       Rcpp::as<Rcpp::DataFrame>(R.parseEval("mysetup$diffusion$vecinj_index"));
// }

void poet::ChemistryParams::initFromR(RInsidePOET &R) {
  // this->database_path =
  //     Rcpp::as<std::string>(R.parseEval("mysetup$chemistry$database"));
  // this->input_script =
  //     Rcpp::as<std::string>(R.parseEval("mysetup$chemistry$input_script"));

  if (R.checkIfExists("dht_species", "mysetup$chemistry")) {
    this->dht_signifs = Rcpp::as<NamedVector<std::uint32_t>>(
        R.parseEval(("mysetup$chemistry$dht_species")));
  }

  if (R.checkIfExists("pht_species", "mysetup$chemistry")) {
    this->pht_signifs = Rcpp::as<NamedVector<std::uint32_t>>(
        R.parseEval(("mysetup$chemistry$pht_species")));
  }
  this->hooks.dht_fill =
      RHookFunction<bool>(R, "mysetup$chemistry$hooks$dht_fill");
  this->hooks.dht_fuzz =
      RHookFunction<std::vector<double>>(R, "mysetup$chemistry$hooks$dht_fuzz");
  this->hooks.interp_pre = RHookFunction<std::vector<std::size_t>>(
      R, "mysetup$chemistry$hooks$interp_pre_func");
  this->hooks.interp_post =
      RHookFunction<bool>(R, "mysetup$chemistry$hooks$interp_post_func");
}

void RuntimeParameters::initVectorParams(RInside &R) {}

std::list<std::string> RuntimeParameters::validateOptions(argh::parser cmdl) {
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
