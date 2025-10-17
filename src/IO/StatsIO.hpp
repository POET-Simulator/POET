#include <string>
#include "Chemistry/ChemistryModule.hpp"

namespace poet
{
    void writeStatsToCSV(const std::vector<ChemistryModule::error_stats> &all_stats,
                         const std::vector<std::string> &species_names,
                         const std::string &out_dir,
                         const std::string &filename);
} // namespace poet
