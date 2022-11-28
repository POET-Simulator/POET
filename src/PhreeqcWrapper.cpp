#include "poet/PhreeqcWrapper.hpp"
#include "IPhreeqc.hpp"
#include "IrmResult.h"
#include "PhreeqcRM.h"

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

poet::PhreeqcWrapper::PhreeqcWrapper(uint32_t inxyz)
    : PhreeqcRM(inxyz, 1), iWPSize(inxyz) {
  this->vecDefMapping.resize(inxyz);

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
  this->RunFile(true, true, false, strInputFile);
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
  // double time = 0.0;
  // double time_step = 0.0;
  // this->SetTime(time);
  // this->SetTimeStep(time_step);
  // this->RunCells();

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

  auto itStart = vecWP.begin();
  auto itEnd = itStart;

  // this->SetMappingForWP(iCurrWPSize);

  int nchem = this->GetChemistryCellCount();

  iCurrElements = this->vecSpeciesPerModule[0];

  itEnd += iCurrElements * this->iWPSize;
  std::vector<double> out;
  this->GetConcentrations(out);
  this->SetConcentrations(std::vector<double>(itStart, itEnd));
  itStart = itEnd;

  // Equlibirum Phases
  if ((iCurrElements = this->vecSpeciesPerModule[1]) != 0) {
    itEnd += iCurrElements * this->iWPSize;
    this->GetPPhaseMoles(out);
    this->SetPPhaseMoles(std::vector<double>(itStart, itEnd));
    itStart = itEnd;
  }

  // // NOTE: Block for 'Surface' and 'Exchange' is missing because of missing
  // // setters @ PhreeqcRM
  // // ...
  // // BLOCK_END

  if ((iCurrElements = this->vecSpeciesPerModule[4]) != 0) {
    itEnd += iCurrElements * this->iWPSize;
    this->GetKineticsMoles(out);
    this->SetKineticsMoles(std::vector<double>(itStart, itEnd));
    itStart = itEnd;
  }
}
void poet::PhreeqcWrapper::GetWPFromInternals(std::vector<double> &vecWP,
                                              uint32_t iCurrWPSize) {
  std::vector<double> vecCurrOutput;

  vecWP.clear();
  vecWP.reserve(this->iSpeciesCount * this->iWPSize);

  this->GetConcentrations(vecCurrOutput);
  vecWP.insert(vecWP.end(), vecCurrOutput.begin(), vecCurrOutput.end());

  if (this->vecSpeciesPerModule[1] != 0) {
    this->GetPPhaseMoles(vecCurrOutput);
    vecWP.insert(vecWP.end(), vecCurrOutput.begin(), vecCurrOutput.end());
  }

  // NOTE: Block for 'Surface' and 'Exchange' is missing because of missing
  // Getters @ PhreeqcRM
  // ...
  // BLOCK_END

  if (this->vecSpeciesPerModule[4] != 0) {
    this->GetKineticsMoles(vecCurrOutput);
    vecWP.insert(vecWP.end(), vecCurrOutput.begin(), vecCurrOutput.end());
  }
}

auto poet::PhreeqcWrapper::ReplaceTotalsByPotentials(
    const std::vector<double> &vecWP, uint32_t iCurrWPSize)
    -> std::vector<double> {

  uint32_t iPropsPerCell = this->vecSpeciesNames.size();

  int iphreeqcResult;
  std::vector<double> vecOutput;
  vecOutput.reserve((iPropsPerCell - 1) * iCurrWPSize);

  auto itCStart = vecWP.begin();
  auto itCEnd = itCStart + (this->vecSpeciesPerModule[0]);
  uint32_t iDiff = iPropsPerCell - this->vecSpeciesPerModule[0];

  for (uint32_t i = 0; i < iCurrWPSize; i++) {
    std::vector<double> vecCIn(itCStart, itCEnd);
    double pH, pe;

    // FIXME: Hardcoded temperatures and pressures here!
    IPhreeqc *util_ptr = this->Concentrations2Utility(
        vecCIn, std::vector<double>(1, 25.0), std::vector<double>(1, 1.0));

    std::string sInput = "SELECTED_OUTPUT " + std::to_string(this->DHT_SELOUT) +
                         "; -pH; -pe;RUN_CELLS; -cells 1; -time_step 0";
    iphreeqcResult = util_ptr->RunString(sInput.c_str());
    if (iphreeqcResult != 0) {
      throw std::runtime_error("IPhreeqc Utility returned non-zero value: " +
                               std::to_string(iphreeqcResult));
    }

    util_ptr->SetCurrentSelectedOutputUserNumber(this->DHT_SELOUT);

    int vtype;
    static std::string svalue(100, '\0');

    iphreeqcResult = util_ptr->GetSelectedOutputValue2(
        1, 0, &vtype, &pH, svalue.data(), svalue.capacity());
    if (iphreeqcResult != 0) {
      throw std::runtime_error("IPhreeqc Utility returned non-zero value: " +
                               std::to_string(iphreeqcResult));
    }

    iphreeqcResult = util_ptr->GetSelectedOutputValue2(
        1, 1, &vtype, &pe, svalue.data(), svalue.capacity());
    if (iphreeqcResult != 0) {
      throw std::runtime_error("IPhreeqc Utility returned non-zero value: " +
                               std::to_string(iphreeqcResult));
    }

    vecOutput.insert(vecOutput.end(), vecCIn.begin() + 3, vecCIn.end());
    vecOutput.push_back(pH);
    vecOutput.push_back(pe);
    vecOutput.insert(vecOutput.end(), itCEnd, itCEnd + iDiff);

    itCStart = itCEnd + iDiff;
    itCEnd = itCStart + (this->vecSpeciesPerModule[0]);
  }

  return vecOutput;
}

IRM_RESULT
poet::PhreeqcWrapper::RunWorkPackage(std::vector<double> &vecWP,
                                     std::vector<int32_t> &vecMapping,
                                     double dSimTime, double dTimestep) {
  if (this->iWPSize != vecMapping.size()) {
    return IRM_INVALIDARG;
  }

  if ((this->iWPSize * this->iSpeciesCount) != vecWP.size()) {
    return IRM_INVALIDARG;
  }

  // check if we actually need to start phreeqc
  bool bRunPhreeqc = false;
  for (const auto &aMappingNum : vecMapping) {
    if (aMappingNum != -1) {
      bRunPhreeqc = true;
      break;
    }
  }

  if (!bRunPhreeqc) {
    return IRM_OK;
  }

  std::vector<double> vecCopy;

  vecCopy = vecWP;
  for (uint32_t i = 0; i < this->iSpeciesCount; i++) {
    for (uint32_t j = 0; j < this->iWPSize; j++) {
      vecWP[(i * this->iWPSize) + j] = vecCopy[(j * this->iSpeciesCount) + i];
    }
  }

  IRM_RESULT result;
  this->CreateMapping(vecMapping);
  this->SetInternalsFromWP(vecWP, this->iWPSize);
  this->SetTime(dSimTime);
  this->SetTimeStep(dTimestep);
  result = this->PhreeqcRM::RunCells();
  this->GetWPFromInternals(vecWP, this->iWPSize);

  vecCopy = vecWP;
  for (uint32_t i = 0; i < this->iSpeciesCount; i++) {
    for (uint32_t j = 0; j < this->iWPSize; j++) {
      vecWP[(j * this->iSpeciesCount) + i] = vecCopy[(i * this->iWPSize) + j];
    }
  }

  return result;
}

inline void poet::PhreeqcWrapper::SetMappingForWP(uint32_t iCurrWPSize) {
  if (iCurrWPSize == this->iWPSize) {
    this->CreateMapping(this->vecDefMapping);
    return;
  }

  std::vector<int> vecCurrMapping(this->vecDefMapping);
  for (uint32_t i = iCurrWPSize; i < vecCurrMapping.size(); i++) {
    vecCurrMapping[i] = -1;
  }
  this->CreateMapping(vecCurrMapping);
}
