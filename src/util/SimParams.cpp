/*
** Copyright (C) 2018-2021 Alexander Lindemann, Max Luebke (University of
** Potsdam)
**
** Copyright (C) 2018-2021 Marco De Lucia (GFZ Potsdam)
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

#include "SimParams.h"

#include <Rcpp.h>

#include <iostream>

using namespace poet;
using namespace std;
using namespace Rcpp;

SimParams::SimParams(int world_rank_, int world_size_) {
  this->simparams.world_rank = world_rank_;
  this->simparams.world_size = world_size_;
}

int SimParams::parseFromCmdl(char *argv[], RRuntime &R) {
  // initialize argh object
  argh::parser cmdl(argv);

  // if user asked for help
  if (cmdl[{"help", "h"}]) {
    if (simparams.world_rank == 0) {
      cout << "Todo" << endl
           << "See README.md for further information." << endl;
    }
    return PARSER_HELP;
  }
  // if positional arguments are missing
  else if (!cmdl(2)) {
    if (simparams.world_rank == 0) {
      cerr << "ERROR. Kin needs 2 positional arguments: " << endl
           << "1) the R script defining your simulation and" << endl
           << "2) the directory prefix where to save results and profiling"
           << endl;
    }
    return PARSER_ERROR;
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
    return PARSER_ERROR;
  }

  /*Parse DHT arguments*/
  simparams.dht_enabled = cmdl["dht"];
  // cout << "CPP: DHT is " << ( dht_enabled ? "ON" : "OFF" ) << '\n';

  if (simparams.dht_enabled) {
    cmdl("dht-strategy", 0) >> simparams.dht_strategy;
    // cout << "CPP: DHT strategy is " << dht_strategy << endl;

    cmdl("dht-signif", 5) >> simparams.dht_significant_digits;
    // cout << "CPP: DHT significant digits = " << dht_significant_digits <<
    // endl;

    simparams.dht_log = !(cmdl["dht-nolog"]);
    // cout << "CPP: DHT logarithm before rounding: " << ( dht_logarithm ? "ON"
    // : "OFF" ) << endl;

    cmdl("dht-size", DHT_SIZE_PER_PROCESS) >> simparams.dht_size_per_process;
    // cout << "CPP: DHT size per process (Byte) = " << dht_size_per_process <<
    // endl;

    cmdl("dht-snaps", 0) >> simparams.dht_snaps;

    cmdl("dht-file") >> dht_file;
  }
  /*Parse work package size*/
  cmdl("work-package-size", WORK_PACKAGE_SIZE_DEFAULT) >> simparams.wp_size;

  /*Parse output options*/
  simparams.store_result = !cmdl["ignore-result"];

  if (simparams.world_rank == 0) {
    cout << "CPP: Complete results storage is "
         << (simparams.store_result ? "ON" : "OFF") << endl;
    cout << "CPP: Work Package Size: " << simparams.wp_size << endl;
    cout << "CPP: DHT is " << (simparams.dht_enabled ? "ON" : "OFF") << '\n';

    if (simparams.dht_enabled) {
      cout << "CPP: DHT strategy is " << simparams.dht_strategy << endl;
      cout << "CPP: DHT key default digits (ignored if 'signif_vector' is "
              "defined) = "
           << simparams.dht_significant_digits << endl;
      cout << "CPP: DHT logarithm before rounding: "
           << (simparams.dht_log ? "ON" : "OFF") << endl;
      cout << "CPP: DHT size per process (Byte) = "
           << simparams.dht_size_per_process << endl;
      cout << "CPP: DHT save snapshots is " << simparams.dht_snaps << endl;
      cout << "CPP: DHT load file is " << dht_file << endl;
    }
  }

  cmdl(1) >> filesim;
  cmdl(2) >> out_dir;

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
  R["dht_enabled"] = simparams.dht_enabled;
  // log before rounding?
  R["dht_log"] = simparams.dht_log;

  // eval the init string, ignoring any returns
  R.parseEvalQ("source(filesim)");

  return PARSER_OK;
}

void SimParams::initVectorParams(RRuntime &R, int col_count) {
  if (simparams.dht_enabled) {
    /*Load significance vector from R setup file (or set default)*/
    bool signif_vector_exists = R.parseEval("exists('signif_vector')");
    if (signif_vector_exists) {
      dht_signif_vector = as<std::vector<int>>(R["signif_vector"]);
    } else {
      dht_signif_vector.assign(col_count, simparams.dht_significant_digits);
    }

    /*Load property type vector from R setup file (or set default)*/
    bool prop_type_vector_exists = R.parseEval("exists('prop_type')");
    if (prop_type_vector_exists) {
      dht_prop_type_vector = as<std::vector<string>>(R["prop_type"]);
    } else {
      dht_prop_type_vector.assign(col_count, "act");
    }

    if (simparams.world_rank == 0) {
      // MDL: new output on signif_vector and prop_type
      if (signif_vector_exists) {
        cout << "CPP: using problem-specific rounding digits: " << endl;
        R.parseEval(
            "print(data.frame(prop=prop, type=prop_type, "
            "digits=signif_vector))");
      } else {
        cout << "CPP: using DHT default rounding digits = "
             << simparams.dht_significant_digits << endl;
      }

      // MDL: pass to R the DHT stuff. These variables exist
      // only if dht_enabled is true
      R["dht_final_signif"] = dht_signif_vector;
      R["dht_final_proptype"] = dht_prop_type_vector;
    }
  }
}

void SimParams::setDtDiffer(bool dt_differ) { simparams.dt_differ = dt_differ; }

t_simparams SimParams::getNumParams() { return this->simparams; }

std::vector<int> SimParams::getDHTSignifVector() {
  return this->dht_signif_vector;
}
std::vector<std::string> SimParams::getDHTPropTypeVector() {
  return this->dht_prop_type_vector;
}
std::string SimParams::getDHTFile() { return this->dht_file; }

std::string SimParams::getFilesim() { return this->filesim; }
std::string SimParams::getOutDir() { return this->out_dir; }

std::list<std::string> SimParams::validateOptions(argh::parser cmdl) {
  /* store all unknown parameters here */
  std::list<std::string> retList;

  /* loop over all flags and compare to given flaglist*/
  for (auto &flag : cmdl.flags()) {
    if (!(flaglist.find(flag) != flaglist.end())) retList.push_back(flag);
  }

  /* and loop also over params and compare to given paramlist */
  for (auto &param : cmdl.params()) {
    if (!(paramlist.find(param.first) != paramlist.end()))
      retList.push_back(param.first);
  }

  return retList;
}
