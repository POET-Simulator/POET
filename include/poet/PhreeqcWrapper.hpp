#ifndef PHREEQCWRAPPER_H_
#define PHREEQCWRAPPER_H_

#include "IrmResult.h"
#include "poet/SimParams.hpp"
#include <PhreeqcRM.h>
#include <cstdint>
#include <string>
#include <vector>

namespace poet {
class PhreeqcWrapper : public PhreeqcRM {
public:
  PhreeqcWrapper(uint32_t inxyz);

  void SetupAndLoadDB(const poet::ChemistryParams &chemPar);
  void InitFromFile(const std::string &strInputFile);

  auto GetSpeciesCount() -> uint32_t;
  auto GetSpeciesCountByModule() -> const std::vector<uint32_t> &;
  auto GetSpeciesNamesFull() -> const std::vector<std::string> &;

  void SetInternalsFromWP(const std::vector<double> &vecWP,
                          uint32_t iCurrWPSize);
  void GetWPFromInternals(std::vector<double> &vecWP, uint32_t iCurrWPSize);

  auto ReplaceTotalsByPotentials(const std::vector<double> &vecWP,
                                 uint32_t iCurrWPSize) -> std::vector<double>;

  IRM_RESULT RunWorkPackage(std::vector<double> &vecWP,
                            std::vector<int32_t> &vecMapping, double dSimTime,
                            double dTimestep);

private:
  void InitInternals();
  inline void SetMappingForWP(uint32_t iCurrWPSize);

  static constexpr int MODULE_COUNT = 5;
  static constexpr uint32_t DHT_SELOUT = 100;

  uint32_t iWPSize;
  bool bInactiveCells = false;
  uint32_t iSpeciesCount;
  std::vector<uint32_t> vecSpeciesPerModule;
  std::vector<std::string> vecSpeciesNames;
  std::vector<int> vecDefMapping;
};
} // namespace poet

#endif // PHREEQCWRAPPER_H_
