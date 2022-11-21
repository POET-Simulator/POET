#include <algorithm>
#include <bits/stdint-uintn.h>
#include <iterator>
#include <poet/PhreeqcWrapper.hpp>
#include <string>
#include <vector>

poet::PhreeqcWrapper::PhreeqcWrapper(uint32_t inxyz)
    : PhreeqcRM(inxyz, 1), iWPSize(inxyz) {
  this->vecDefMapping.reserve(inxyz);

  for (uint32_t i = 0; i < inxyz; i++) {
    this->vecDefMapping[i] = i;
  }
}

void poet::PhreeqcWrapper::SetupAndLoadDB(
    const poet::ChemistryParams &chemPar) {
  // TODO: hardcoded options ...
  this->SetErrorHandlerMode(1);
  this->SetComponentH2O(false);
  this->SetRebalanceFraction(0.5);
  this->SetRebalanceByCell(true);
  this->UseSolutionDensityVolume(false);
  this->SetPartitionUZSolids(false);

  // Set concentration units
  // 1, mg/L; 2, mol/L; 3, kg/kgs
  this->SetUnitsSolution(2);
  // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
  this->SetUnitsPPassemblage(1);
  // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
  this->SetUnitsExchange(1);
  // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
  this->SetUnitsSurface(1);
  // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
  this->SetUnitsGasPhase(1);
  // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
  this->SetUnitsSSassemblage(1);
  // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
  this->SetUnitsKinetics(1);

  // Set representative volume
  std::vector<double> rv;
  rv.resize(this->iWPSize, 1.0);
  this->SetRepresentativeVolume(rv);

  // Set initial porosity
  std::vector<double> por;
  por.resize(this->iWPSize, 0.05);
  this->SetPorosity(por);

  // Set initial saturation
  std::vector<double> sat;
  sat.resize(this->iWPSize, 1.0);
  this->SetSaturation(sat);

  // Load database
  this->LoadDatabase(chemPar.database_path);
}

void poet::PhreeqcWrapper::InitFromFile(const std::string &strInputFile) {
  this->RunFile(true, true, true, strInputFile);
  this->RunString(true, false, true, "DELETE; -all");

  this->FindComponents();

  std::vector<int> ic1;
  ic1.resize(this->iWPSize * 7, -1);
  // TODO: hardcoded reaction modules
  for (int i = 0; i < nxyz; i++) {
    ic1[i] = 1;             // Solution 1
    ic1[nxyz + i] = 1;      // Equilibrium 1
    ic1[2 * nxyz + i] = -1; // Exchange none
    ic1[3 * nxyz + i] = -1; // Surface none
    ic1[4 * nxyz + i] = -1; // Gas phase none
    ic1[5 * nxyz + i] = -1; // Solid solutions none
    ic1[6 * nxyz + i] = 1;  // Kinetics 1
  }

  this->InitialPhreeqc2Module(ic1);

  // Initial equilibration of cells
  double time = 0.0;
  double time_step = 0.0;
  this->SetTime(time);
  this->SetTimeStep(time_step);
  this->RunCells();

  this->InitInternals();
}

void poet::PhreeqcWrapper::InitInternals() {
  uint32_t sum = 0, curr_count;
  std::vector<std::string> vecCurrSpecies;

  this->vecSpeciesPerModule.reserve(this->MODULE_COUNT);

  vecCurrSpecies = this->GetComponents();
  this->vecSpeciesNames.insert(this->vecSpeciesNames.end(),
                               vecCurrSpecies.begin(), vecCurrSpecies.end());
  this->vecSpeciesPerModule.push_back(vecCurrSpecies.size());

  vecCurrSpecies = this->GetEquilibriumPhases();
  this->vecSpeciesNames.insert(this->vecSpeciesNames.end(),
                               vecCurrSpecies.begin(), vecCurrSpecies.end());
  this->vecSpeciesPerModule.push_back(vecCurrSpecies.size());

  vecCurrSpecies = this->GetExchangeSpecies();
  this->vecSpeciesNames.insert(this->vecSpeciesNames.end(),
                               vecCurrSpecies.begin(), vecCurrSpecies.end());
  this->vecSpeciesPerModule.push_back(vecCurrSpecies.size());

  vecCurrSpecies = this->GetSurfaceSpecies();
  this->vecSpeciesNames.insert(this->vecSpeciesNames.end(),
                               vecCurrSpecies.begin(), vecCurrSpecies.end());
  this->vecSpeciesPerModule.push_back(vecCurrSpecies.size());

  vecCurrSpecies = this->GetKineticReactions();
  this->vecSpeciesNames.insert(this->vecSpeciesNames.end(),
                               vecCurrSpecies.begin(), vecCurrSpecies.end());
  this->vecSpeciesPerModule.push_back(vecCurrSpecies.size());

  this->iSpeciesCount = this->vecSpeciesNames.size();
}

auto poet::PhreeqcWrapper::GetSpeciesCount() -> uint32_t {
  return this->iSpeciesCount;
}

auto poet::PhreeqcWrapper::GetSpeciesCountByModule()
    -> const std::vector<uint32_t> & {
  return this->vecSpeciesPerModule;
}

auto poet::PhreeqcWrapper::GetSpeciesNamesFull()
    -> const std::vector<std::string> & {
  return this->vecSpeciesNames;
}

void poet::PhreeqcWrapper::SetInternalsFromWP(const std::vector<double> &vecWP,
                                              uint32_t iCurrWPSize) {
  uint32_t iCurrElements;

  if (iCurrWPSize != this->iWPSize) {
    this->SetMappingForWP(iCurrWPSize);
  }

  auto itStart = vecWP.begin();
  auto itEnd = itStart;

  iCurrElements = this->vecSpeciesPerModule[0];

  itEnd += iCurrElements * iCurrWPSize;
  this->SetConcentrations(std::vector<double>(itStart, itEnd));
  itStart = itEnd;

  // Equlibirum Phases
  if ((iCurrElements = this->vecSpeciesPerModule[1]) != 0) {
    itEnd += iCurrElements * iCurrWPSize;
    this->SetPPhaseMoles(std::vector<double>(itStart, itEnd));
    itStart = itEnd;
  }

  // // NOTE: Block for 'Surface' and 'Exchange' is missing because of missing
  // // setters @ PhreeqcRM
  // // ...
  // // BLOCK_END

  if ((iCurrElements = this->vecSpeciesPerModule[4]) != 0) {
    itEnd += iCurrElements * iCurrWPSize;
    this->SetKineticsMoles(std::vector<double>(itStart, itEnd));
    itStart = itEnd;
  }
}
auto poet::PhreeqcWrapper::GetWPFromInternals(uint32_t iCurrWPSize)
    -> std::vector<double> {
  std::vector<double> vecOutput, vecCurrOutput;

  vecOutput.reserve(this->iSpeciesCount * iCurrWPSize);

  this->GetConcentrations(vecCurrOutput);
  vecOutput.insert(vecOutput.end(), vecCurrOutput.begin(), vecCurrOutput.end());

  if (this->vecSpeciesPerModule[1] != 0) {
    this->GetPPhaseMoles(vecCurrOutput);
    vecOutput.insert(vecOutput.end(), vecCurrOutput.begin(),
                     vecCurrOutput.end());
  }

  // NOTE: Block for 'Surface' and 'Exchange' is missing because of missing
  // Getters @ PhreeqcRM
  // ...
  // BLOCK_END

  if (this->vecSpeciesPerModule[4] != 0) {
    this->GetKineticsMoles(vecCurrOutput);
    vecOutput.insert(vecOutput.end(), vecCurrOutput.begin(),
                     vecCurrOutput.end());
  }

  return vecOutput;
}

void poet::PhreeqcWrapper::SetMappingForWP(uint32_t iCurrWPSize) {
  if (iCurrWPSize == this->iWPSize) {
    this->CreateMapping(this->vecDefMapping);
  } else {
    std::vector<int> vecCurrMapping = this->vecDefMapping;
    for (uint32_t i = iCurrWPSize; i < vecCurrMapping.size(); i++) {
      vecCurrMapping[i] = -1;
    }
    this->CreateMapping(vecCurrMapping);
  }
}
