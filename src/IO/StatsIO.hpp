#include <string>
#include "Chemistry/ChemistryModule.hpp"

namespace poet
{
    void writeStatsToCSV(const std::vector<ChemistryModule::SimulationErrorStats> &all_stats,
                         const std::vector<std::string> &species_names,
                         const std::string &filename);
} // namespace poet
