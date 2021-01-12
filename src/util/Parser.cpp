#include "Parser.h"

#include <Rcpp.h>

#include <iostream>

using namespace poet;
using namespace std;
using namespace Rcpp;

Parser::Parser(char *argv[], int world_rank_, int world_size_)
    : cmdl(argv), world_rank(world_rank_), world_size(world_size_) {
  this->simparams.world_rank = world_rank_;
  this->simparams.world_size = world_size;
}

int Parser::parseCmdl() {
  // if user asked for help
  if (cmdl[{"help", "h"}]) {
    if (world_rank == 0) {
      cout << "Todo" << endl
           << "See README.md for further information." << endl;
    }
    return PARSER_HELP;
  }
  // if positional arguments are missing
  else if (!cmdl(2)) {
    if (world_rank == 0) {
      cerr << "ERROR. Kin needs 2 positional arguments: " << endl
           << "1) the R script defining your simulation and" << endl
           << "2) the directory prefix where to save results and profiling"
           << endl;
    }
    return PARSER_ERROR;
  }

  std::list<std::string> optionsError = checkOptions(cmdl);
  if (!optionsError.empty()) {
    if (world_rank == 0) {
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

    cmdl("dht-signif", 5) >> dht_significant_digits;
    // cout << "CPP: DHT significant digits = " << dht_significant_digits <<
    // endl;

    simparams.dht_log = !(cmdl["dht-nolog"]);
    // cout << "CPP: DHT logarithm before rounding: " << ( dht_logarithm ? "ON"
    // : "OFF" ) << endl;

    cmdl("dht-size", DHT_SIZE_PER_PROCESS) >> simparams.dht_size_per_process;
    // cout << "CPP: DHT size per process (Byte) = " << dht_size_per_process <<
    // endl;

    cmdl("dht-snaps", 0) >> simparams.dht_snaps;

    cmdl("dht-file") >> simparams.dht_file;
  }
  /*Parse work package size*/
  cmdl("work-package-size", WORK_PACKAGE_SIZE_DEFAULT) >> simparams.wp_size;

  /*Parse output options*/
  simparams.store_result = !cmdl["ignore-result"];

  if (world_rank == 0) {
    cout << "CPP: Complete results storage is "
         << (simparams.store_result ? "ON" : "OFF") << endl;
    cout << "CPP: Work Package Size: " << simparams.wp_size << endl;
    cout << "CPP: DHT is " << (simparams.dht_enabled ? "ON" : "OFF") << '\n';

    if (simparams.dht_enabled) {
      cout << "CPP: DHT strategy is " << simparams.dht_strategy << endl;
      cout << "CPP: DHT key default digits (ignored if 'signif_vector' is "
              "defined) = "
           << dht_significant_digits << endl;
      cout << "CPP: DHT logarithm before rounding: "
           << (simparams.dht_log ? "ON" : "OFF") << endl;
      cout << "CPP: DHT size per process (Byte) = "
           << simparams.dht_size_per_process << endl;
      cout << "CPP: DHT save snapshots is " << simparams.dht_snaps << endl;
      cout << "CPP: DHT load file is " << simparams.dht_file << endl;
    }
  }

  cmdl(1) >> simparams.filesim;
  cmdl(2) >> simparams.out_dir;

  return PARSER_OK;
}

void Parser::parseR(RRuntime &R) {
  // if local_rank == 0 then master else worker
  R["local_rank"] = simparams.world_rank;
  // assign a char* (string) to 'filesim'
  R["filesim"] = wrap(simparams.filesim);
  // assign a char* (string) to 'fileout'
  R["fileout"] = wrap(simparams.out_dir);
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
}

t_simparams Parser::getParams() { return this->simparams; }

std::list<std::string> Parser::checkOptions(argh::parser cmdl) {
  std::list<std::string> retList;
  //   std::set<std::string> flist = flagList();
  //   std::set<std::string> plist = paramList();

  for (auto &flag : cmdl.flags()) {
    if (!(flaglist.find(flag) != flaglist.end())) retList.push_back(flag);
  }

  for (auto &param : cmdl.params()) {
    if (!(paramlist.find(param.first) != paramlist.end()))
      retList.push_back(param.first);
  }

  return retList;
}
