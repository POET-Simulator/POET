#include "ChemistryModule.hpp"

#include "PhreeqcEngine.hpp"
#include "PhreeqcMatrix.hpp"
#include "SurrogateModels/DHT_Wrapper.hpp"
#include "SurrogateModels/Interpolation.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <memory>
#include <mpi.h>
#include <string>
#include <vector>

std::vector<double>
inverseDistanceWeighting(const std::vector<std::int32_t> &to_calc,
                         const std::vector<double> &from,
                         const std::vector<std::vector<double>> &input,
                         const std::vector<std::vector<double>> &output) {
  std::vector<double> results = from;

  // const std::uint32_t buffer_size = input.size() + 1;
  // double buffer[buffer_size];
  // double from_rescaled;

  const std::uint32_t data_set_n = input.size();
  double rescaled[to_calc.size()][data_set_n + 1];
  double weights[data_set_n];

  // rescaling over all key elements
  for (std::uint32_t key_comp_i = 0; key_comp_i < to_calc.size();
       key_comp_i++) {
    const auto output_comp_i = to_calc[key_comp_i];

    // rescale input between 0 and 1
    for (std::uint32_t point_i = 0; point_i < data_set_n; point_i++) {
      rescaled[key_comp_i][point_i] = input[point_i][key_comp_i];
    }

    rescaled[key_comp_i][data_set_n] = from[output_comp_i];

    const double min = *std::min_element(rescaled[key_comp_i],
                                         rescaled[key_comp_i] + data_set_n + 1);
    const double max = *std::max_element(rescaled[key_comp_i],
                                         rescaled[key_comp_i] + data_set_n + 1);

    for (std::uint32_t point_i = 0; point_i < data_set_n; point_i++) {
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
  for (std::uint32_t point_i = 0; point_i < data_set_n; point_i++) {
    double distance = 0;
    for (std::uint32_t key_comp_i = 0; key_comp_i < to_calc.size();
         key_comp_i++) {
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

  for (std::uint32_t key_comp_i = 0; key_comp_i < to_calc.size();
       key_comp_i++) {
    const auto output_comp_i = to_calc[key_comp_i];
    double key_delta = 0;

    // if (interp_i == 0) {
    //   has_h = true;
    // }

    // if (interp_i == 1) {
    //   has_o = true;
    // }

    for (std::uint32_t j = 0; j < data_set_n; j++) {
      key_delta += weights[j] * output[j][output_comp_i];
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

poet::ChemistryModule::ChemistryModule(
    uint32_t wp_size_, const InitialList::ChemistryInit chem_params,
    MPI_Comm communicator)
    : group_comm(communicator), wp_size(wp_size_), params(chem_params) {
  MPI_Comm_rank(communicator, &comm_rank);
  MPI_Comm_size(communicator, &comm_size);

  this->is_sequential = comm_size == 1;
  this->is_master = comm_rank == 0;

  this->n_cells = chem_params.total_grid_cells;

  if (!is_master) {
    PhreeqcMatrix pqc_mat =
        PhreeqcMatrix(chem_params.database, chem_params.pqc_script);
    for (const auto &cell_id : chem_params.pqc_ids) {
      this->phreeqc_instances[cell_id] =
          std::make_unique<PhreeqcEngine>(pqc_mat, cell_id);
    }
    // for (std::size_t i = 0; i < chem_params.pqc_ids.size(); i++) {
    //   this->phreeqc_instances[chem_params.pqc_ids[i]] =
    //       std::make_unique<PhreeqcWrapper>(
    //           chem_params.database, chem_params.pqc_scripts[i],
    //           chem_params.pqc_sol_order, chem_params.field_header, wp_size_);
    // }
  }
}

poet::ChemistryModule::~ChemistryModule() {
  if (dht) {
    delete dht;
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
      std::vector<std::uint32_t> default_signif(
          this->prop_count, DHT_Wrapper::DHT_KEY_SIGNIF_DEFAULT);
      map_copy = NamedVector<std::uint32_t>(this->prop_names, default_signif);
    }

    auto key_indices = parseDHTSpeciesVec(key_species, this->prop_names);

    if (this->dht) {
      delete this->dht;
    }

    const std::uint64_t dht_size = size_mb * MB_FACTOR;

    this->dht = new DHT_Wrapper(dht_comm, dht_size, map_copy, key_indices,
                                this->prop_names, params.hooks,
                                this->prop_count, interp_enabled);
    this->dht->setBaseTotals(base_totals.at(0), base_totals.at(1));
  }
}

inline std::vector<std::int32_t> poet::ChemistryModule::parseDHTSpeciesVec(
    const NamedVector<std::uint32_t> &key_species,
    const std::vector<std::string> &to_compare) const {
  std::vector<int32_t> species_indices;
  species_indices.reserve(key_species.size());

  const auto test = key_species.getNames();

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
      for (auto i = 0; i < map_copy.size(); i++) {
        const std::uint32_t signif =
            static_cast<std::uint32_t>(map_copy[i]) - (map_copy[i] > InterpolationModule::COARSE_DIFF
                               ? InterpolationModule::COARSE_DIFF
                               : 0);
        map_copy[i] = signif;
      }
    }

    auto key_indices =
        parseDHTSpeciesVec(map_copy, dht->getKeySpecies().getNames());

    if (this->interp) {
      this->interp.reset();
    }

    const uint64_t pht_size = size_mb * MB_FACTOR;

    interp = std::make_unique<poet::InterpolationModule>(
        bucket_size, pht_size, min_entries, *(this->dht), map_copy, key_indices,
        this->prop_names, this->params.hooks);

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
  
void poet::ChemistryModule::set_ai_surrogate_validity_vector(std::vector<int> r_vector) {
  this->ai_surrogate_validity_vector = r_vector;
}
