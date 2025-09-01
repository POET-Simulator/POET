#include "IO/StatsIO.hpp"
#include <fstream>
#include <iostream>
#include <string>

namespace poet
{
    void writeStatsToCSV(const std::vector<ChemistryModule::error_stats> &all_stats,
                         const std::vector<std::string> &species_names,
                         const std::string &filename)
    {
        std::ofstream out(filename);
        if (!out.is_open())
        {
            std::cerr << "Could not open " << filename << " !" << std::endl;
            return;
        }

        // header
        out << "Iteration, Species, MAPE, RRSME \n";

        for (size_t i = 0; i < all_stats.size(); ++i)
        {
            for (size_t j = 0; j < species_names.size(); ++j)
            {
                out << all_stats[i].iteration << ",\t"
                    << species_names[j] << ",\t"
                    << all_stats[i].mape[j] << ",\t"
                    << all_stats[i].rrsme[j] << "\n";
            }
            out << std::endl;
        }

        out.close();
        std::cout << "Stats written to " << filename << "\n";
    }
} // namespace poet