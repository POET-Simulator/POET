#include "poet/ChemistryModule.hpp"
#include "PhreeqcRM.h"
#include "poet/DHT_Wrapper.hpp"
#include "poet/Interpolation.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <mpi.h>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

constexpr uint32_t MB_FACTOR = 1E6;

std::vector<double>
inverseDistanceWeighting(const std::vector<std::int32_t> &to_calc,
                         const std::vector<double> &from,
                         const std::vector<std::vector<double>> &input,
                         const std::vector<std::vector<double>> &output) {
  std::vector<double> results = from;

  const std::uint32_t buffer_size = input.size() + 1;
  double buffer[buffer_size];
  double from_rescaled;

  const std::uint32_t data_set_n = input.size();
  double rescaled[to_calc.size()][data_set_n + 1];
  double weights[data_set_n];

  // rescaling over all key elements
  for (int key_comp_i = 0; key_comp_i < to_calc.size(); key_comp_i++) {
    const auto output_comp_i = to_calc[key_comp_i];

    // rescale input between 0 and 1
    for (int point_i = 0; point_i < data_set_n; point_i++) {
      rescaled[key_comp_i][point_i] = input[point_i][key_comp_i];
    }

    rescaled[key_comp_i][data_set_n] = from[output_comp_i];

    const double min = *std::min_element(rescaled[key_comp_i],
                                         rescaled[key_comp_i] + data_set_n + 1);
    const double max = *std::max_element(rescaled[key_comp_i],
                                         rescaled[key_comp_i] + data_set_n + 1);

    for (int point_i = 0; point_i < data_set_n; point_i++) {
      rescaled[key_comp_i][point_i] =
          ((max - min) != 0
               ? (rescaled[key_comp_i][point_i] - min) / (max - min)
               : 0);
    }
    rescaled[key_comp_i][data_set_n] =
        ((max - min) != 0 ? (from[output_comp_i] - min) / (max - min) : 0);
  }

  // calculate distances for each data set
  double inv_sum = 0;
  for (int point_i = 0; point_i < data_set_n; point_i++) {
    double distance = 0;
    for (int key_comp_i = 0; key_comp_i < to_calc.size(); key_comp_i++) {
      distance += std::pow(
          rescaled[key_comp_i][point_i] - rescaled[key_comp_i][data_set_n], 2);
    }
    weights[point_i] = 1 / std::sqrt(distance);
    assert(!std::isnan(weights[point_i]));
    inv_sum += weights[point_i];
  }

  assert(!std::isnan(inv_sum));

  // actual interpolation
  // bool has_h = false;
  // bool has_o = false;

  for (int key_comp_i = 0; key_comp_i < to_calc.size(); key_comp_i++) {
    const auto output_comp_i = to_calc[key_comp_i];
    double key_delta = 0;

    // if (interp_i == 0) {
    //   has_h = true;
    // }

    // if (interp_i == 1) {
    //   has_o = true;
    // }

    for (int j = 0; j < data_set_n; j++) {
      key_delta += weights[j] * output[j][key_comp_i];
    }

    key_delta /= inv_sum;

    results[output_comp_i] = from[output_comp_i] + key_delta;
  }

  // if (!has_h) {
  //   double new_val = 0;
  //   for (int j = 0; j < data_set_n; j++) {
  //     new_val += weights[j] * output[j][0];
  //   }
  //   results[0] = new_val / inv_sum;
  // }

  // if (!has_h) {
  //   double new_val = 0;
  //   for (int j = 0; j < data_set_n; j++) {
  //     new_val += weights[j] * output[j][1];
  //   }
  //   results[1] = new_val / inv_sum;
  // }

  // for (std::uint32_t i = 0; i < to_calc.size(); i++) {
  //   const std::uint32_t interp_i = to_calc[i];

  //   // rescale input between 0 and 1
  //   for (int j = 0; j < input.size(); j++) {
  //     buffer[j] = input[j].at(i);
  //   }

  //   buffer[buffer_size - 1] = from[interp_i];

  //   const double min = *std::min_element(buffer, buffer + buffer_size);
  //   const double max = *std::max_element(buffer, buffer + buffer_size);

  //   for (int j = 0; j < input.size(); j++) {
  //     buffer[j] = ((max - min) != 0 ? (buffer[j] - min) / (max - min) : 1);
  //   }
  //   from_rescaled =
  //       ((max - min) != 0 ? (from[interp_i] - min) / (max - min) : 0);

  //   double inv_sum = 0;

  //   // calculate distances for each point
  //   for (int i = 0; i < input.size(); i++) {
  //     const double distance = std::pow(buffer[i] - from_rescaled, 2);

  //     buffer[i] = distance > 0 ? (1 / std::sqrt(distance)) : 0;
  //     inv_sum += buffer[i];
  //   }
  //   // calculate new values
  //   double new_val = 0;
  //   for (int i = 0; i < output.size(); i++) {
  //     new_val += buffer[i] * output[i][interp_i];
  //   }
  //   results[interp_i] = new_val / inv_sum;
  //   if (std::isnan(results[interp_i])) {
  //     std::cout << "nan with new_val = " << output[0][i] << std::endl;
  //   }
  // }

  return results;
}

poet::ChemistryModule::ChemistryModule(uint32_t nxyz, uint32_t wp_size,
                                       std::uint32_t maxiter,
                                       const ChemistryParams &chem_param,
                                       MPI_Comm communicator)
    : PhreeqcRM(nxyz, 1), group_comm(communicator), wp_size(wp_size),
      params(chem_param) {

  MPI_Comm_size(communicator, &this->comm_size);
  MPI_Comm_rank(communicator, &this->comm_rank);

  this->is_sequential = (this->comm_size == 1);
  this->is_master = (this->comm_rank == 0);

  if (!is_sequential && is_master) {
    MPI_Bcast(&wp_size, 1, MPI_UINT32_T, 0, this->group_comm);
    MPI_Bcast(&maxiter, 1, MPI_UINT32_T, 0, this->group_comm);
  }

  this->file_pad = std::ceil(std::log10(maxiter + 1));
}

poet::ChemistryModule::~ChemistryModule() {
  if (dht) {
    delete dht;
  }
}

poet::ChemistryModule
poet::ChemistryModule::createWorker(MPI_Comm communicator,
                                    const ChemistryParams &chem_param) {
  uint32_t wp_size;
  MPI_Bcast(&wp_size, 1, MPI_UINT32_T, 0, communicator);

  std::uint32_t maxiter;
  MPI_Bcast(&maxiter, 1, MPI_UINT32_T, 0, communicator);

  return ChemistryModule(wp_size, wp_size, maxiter, chem_param, communicator);
}

void poet::ChemistryModule::RunInitFile(const std::string &input_script_path) {
  if (this->is_master) {
    int f_type = CHEM_INIT;
    PropagateFunctionType(f_type);

    int count = input_script_path.size();
    ChemBCast(&count, 1, MPI_INT);
    ChemBCast(const_cast<char *>(input_script_path.data()), count, MPI_CHAR);
  }

  this->RunFile(true, true, false, input_script_path);
  this->RunString(true, false, false, "DELETE; -all; PRINT; -warnings 0;");

  this->FindComponents();

  PhreeqcRM::initializePOET(this->speciesPerModule, this->prop_names);
  this->prop_count = prop_names.size();

  char exchange = (speciesPerModule[1] == 0 ? -1 : 1);
  char kinetics = (speciesPerModule[2] == 0 ? -1 : 1);
  char equilibrium = (speciesPerModule[3] == 0 ? -1 : 1);
  char surface = (speciesPerModule[4] == 0 ? -1 : 1);

  std::vector<int> ic1;
  ic1.resize(this->nxyz * 7, -1);
  // TODO: hardcoded reaction modules
  for (int i = 0; i < nxyz; i++) {
    ic1[i] = 1;                   // Solution 1
    ic1[nxyz + i] = equilibrium;  // Equilibrium 1
    ic1[2 * nxyz + i] = exchange; // Exchange none
    ic1[3 * nxyz + i] = surface;  // Surface none
    ic1[4 * nxyz + i] = -1;       // Gas phase none
    ic1[5 * nxyz + i] = -1;       // Solid solutions none
    ic1[6 * nxyz + i] = kinetics; // Kinetics 1
  }

  this->InitialPhreeqc2Module(ic1);
}

void poet::ChemistryModule::initializeField(const Field &trans_field) {

  if (is_master) {
    int f_type = CHEM_INIT_SPECIES;
    PropagateFunctionType(f_type);
  }

  std::vector<std::string> essentials_backup{
      prop_names.begin() + speciesPerModule[0], prop_names.end()};

  std::vector<std::string> new_solution_names{
      this->prop_names.begin(), this->prop_names.begin() + speciesPerModule[0]};

  if (is_master) {
    for (auto &prop : trans_field.GetProps()) {
      if (std::find(new_solution_names.begin(), new_solution_names.end(),
                    prop) == new_solution_names.end()) {
        int size = prop.size();
        ChemBCast(&size, 1, MPI_INT);
        ChemBCast(prop.data(), prop.size(), MPI_CHAR);
        new_solution_names.push_back(prop);
      }
    }
    int end = 0;
    ChemBCast(&end, 1, MPI_INT);
  } else {
    constexpr int MAXSIZE = 128;
    MPI_Status status;
    int recv_size;
    char recv_buffer[MAXSIZE];
    while (1) {
      ChemBCast(&recv_size, 1, MPI_INT);
      if (recv_size == 0) {
        break;
      }
      ChemBCast(recv_buffer, recv_size, MPI_CHAR);
      recv_buffer[recv_size] = '\0';
      new_solution_names.push_back(std::string(recv_buffer));
    }
  }

  // now sort the new values
  std::sort(new_solution_names.begin() + 3, new_solution_names.end());
  this->SetPOETSolutionNames(new_solution_names);
  this->speciesPerModule[0] = new_solution_names.size();

  // and append other processes than solutions
  std::vector<std::string> new_prop_names = new_solution_names;
  new_prop_names.insert(new_prop_names.end(), essentials_backup.begin(),
                        essentials_backup.end());

  std::vector<std::string> old_prop_names{this->prop_names};
  this->prop_names = std::move(new_prop_names);
  this->prop_count = prop_names.size();

  if (is_master) {
    this->n_cells = trans_field.GetRequestedVecSize();
    chem_field = Field(n_cells);

    std::vector<double> phreeqc_init;
    this->getDumpedField(phreeqc_init);

    if (is_sequential) {
      std::vector<double> init_vec{phreeqc_init};
      this->unshuffleField(phreeqc_init, n_cells, prop_count, 1, init_vec);
      chem_field.InitFromVec(init_vec, prop_names);
      return;
    }

    std::vector<std::vector<double>> initial_values;

    for (const auto &vec : trans_field.As2DVector()) {
      initial_values.push_back(vec);
    }

    this->base_totals = {initial_values.at(0).at(0),
                         initial_values.at(1).at(0)};
    ChemBCast(base_totals.data(), 2, MPI_DOUBLE);

    for (int i = speciesPerModule[0]; i < phreeqc_init.size(); i++) {
      std::vector<double> init(n_cells, phreeqc_init[i]);
      initial_values.push_back(std::move(init));
    }

    chem_field.InitFromVec(initial_values, prop_names);
  } else {
    ChemBCast(base_totals.data(), 2, MPI_DOUBLE);
  }

  if (this->params.use_dht || this->params.use_interp) {
    initializeDHT(this->params.dht_size, this->params.dht_signifs);
    setDHTSnapshots(this->params.dht_snaps, this->params.dht_outdir);
    setDHTReadFile(this->params.dht_file);

    this->dht_enabled = this->params.use_dht;

    if (this->params.use_interp) {
      initializeInterp(this->params.pht_max_entries, this->params.pht_size,
                       this->params.interp_min_entries,
                       this->params.pht_signifs);
      this->interp_enabled = this->params.use_interp;
    }
  }
}

void poet::ChemistryModule::initializeDHT(
    uint32_t size_mb, const NamedVector<std::uint32_t> &key_species) {
  constexpr uint32_t MB_FACTOR = 1E6;

  this->dht_enabled = true;

  MPI_Comm dht_comm;

  if (this->is_master) {
    MPI_Comm_split(this->group_comm, MPI_UNDEFINED, this->comm_rank, &dht_comm);
    return;
  }

  if (!this->is_master) {

    MPI_Comm_split(this->group_comm, 1, this->comm_rank, &dht_comm);

    auto map_copy = key_species;

    if (key_species.empty()) {
      const auto &all_species = this->prop_names;
      for (std::size_t i = 0; i < all_species.size(); i++) {
        map_copy.insert(std::make_pair(all_species[i],
                                       DHT_Wrapper::DHT_KEY_SIGNIF_DEFAULT));
      }
    }

    auto key_indices = parseDHTSpeciesVec(key_species, this->prop_names);

    if (this->dht) {
      delete this->dht;
    }

    const std::uint64_t dht_size = size_mb * MB_FACTOR;

    this->dht = new DHT_Wrapper(dht_comm, dht_size, map_copy, key_indices,
                                this->prop_count);
    this->dht->setBaseTotals(base_totals.at(0), base_totals.at(1));
  }
}

inline std::vector<std::int32_t> poet::ChemistryModule::parseDHTSpeciesVec(
    const NamedVector<std::uint32_t> &key_species,
    const std::vector<std::string> &to_compare) const {
  std::vector<int32_t> species_indices;
  species_indices.reserve(key_species.size());

  for (const auto &species : key_species.getNames()) {
    auto it = std::find(to_compare.begin(), to_compare.end(), species);
    if (it == to_compare.end()) {
      species_indices.push_back(DHT_Wrapper::DHT_KEY_INPUT_CUSTOM);
      continue;
    }
    const std::uint32_t index = it - to_compare.begin();
    species_indices.push_back(index);
  }

  return species_indices;
}

void poet::ChemistryModule::BCastStringVec(std::vector<std::string> &io) {

  if (this->is_master) {
    int vec_size = io.size();
    ChemBCast(&vec_size, 1, MPI_INT);

    for (const auto &value : io) {
      int buf_size = value.size() + 1;
      ChemBCast(&buf_size, 1, MPI_INT);
      ChemBCast(const_cast<char *>(value.c_str()), buf_size, MPI_CHAR);
    }
  } else {
    int vec_size;
    ChemBCast(&vec_size, 1, MPI_INT);

    io.resize(vec_size);

    for (int i = 0; i < vec_size; i++) {
      int buf_size;
      ChemBCast(&buf_size, 1, MPI_INT);
      char buf[buf_size];
      ChemBCast(buf, buf_size, MPI_CHAR);
      io[i] = std::string{buf};
    }
  }
}

void poet::ChemistryModule::setDHTSnapshots(int type,
                                            const std::string &out_dir) {
  if (this->is_master) {
    return;
  }

  this->dht_file_out_dir = out_dir;
  this->dht_snaps_type = type;
}

void poet::ChemistryModule::setDHTReadFile(const std::string &input_file) {
  if (this->is_master) {
    return;
  }

  if (!input_file.empty()) {
    WorkerReadDHTDump(input_file);
  }
}

void poet::ChemistryModule::initializeInterp(
    std::uint32_t bucket_size, std::uint32_t size_mb, std::uint32_t min_entries,
    const NamedVector<std::uint32_t> &key_species) {

  if (!this->is_master) {

    constexpr uint32_t MB_FACTOR = 1E6;

    assert(this->dht);

    this->interp_enabled = true;

    auto map_copy = key_species;

    if (key_species.empty()) {
      map_copy = this->dht->getKeySpecies();
      for (std::size_t i = 0; i < map_copy.size(); i++) {
        const std::uint32_t signif =
            map_copy[i] - (map_copy[i] > InterpolationModule::COARSE_DIFF
                               ? InterpolationModule::COARSE_DIFF
                               : 0);
        map_copy[i] = signif;
      }
    }

    auto key_indices =
        parseDHTSpeciesVec(key_species, dht->getKeySpecies().getNames());

    if (this->interp) {
      this->interp.reset();
    }

    const uint64_t pht_size = size_mb * MB_FACTOR;

    interp = std::make_unique<poet::InterpolationModule>(
        bucket_size, pht_size, min_entries, *(this->dht), map_copy,
        key_indices);

    interp->setInterpolationFunction(inverseDistanceWeighting);
  }
}

std::vector<double>
poet::ChemistryModule::shuffleField(const std::vector<double> &in_field,
                                    uint32_t size_per_prop, uint32_t prop_count,
                                    uint32_t wp_count) {
  std::vector<double> out_buffer(in_field.size());
  uint32_t write_i = 0;
  for (uint32_t i = 0; i < wp_count; i++) {
    for (uint32_t j = i; j < size_per_prop; j += wp_count) {
      for (uint32_t k = 0; k < prop_count; k++) {
        out_buffer[(write_i * prop_count) + k] =
            in_field[(k * size_per_prop) + j];
      }
      write_i++;
    }
  }
  return out_buffer;
}

void poet::ChemistryModule::unshuffleField(const std::vector<double> &in_buffer,
                                           uint32_t size_per_prop,
                                           uint32_t prop_count,
                                           uint32_t wp_count,
                                           std::vector<double> &out_field) {
  uint32_t read_i = 0;

  for (uint32_t i = 0; i < wp_count; i++) {
    for (uint32_t j = i; j < size_per_prop; j += wp_count) {
      for (uint32_t k = 0; k < prop_count; k++) {
        out_field[(k * size_per_prop) + j] =
            in_buffer[(read_i * prop_count) + k];
      }
      read_i++;
    }
  }
}
