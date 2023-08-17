#include "poet/DataStructures.hpp"
#include <Rcpp.h>
#include <Rinternals.h>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

poet::Field::Field(std::uint32_t vec_s,
                   const std::vector<std::vector<double>> &input,
                   const std::vector<std::string> &prop_vec, bool row_major)
    : req_vec_size(vec_s), props(prop_vec) {
  if (row_major) {
    if (this->props.size() != input[0].size()) {
      throw std::runtime_error("Prop vector shall name all elements.");
    }

    const std::size_t n_input = input.size();

    for (std::size_t i = 0; i < this->props.size(); i++) {
      std::vector<double> curr_col(n_input);

      for (std::size_t j = 0; j < n_input; j++) {
        curr_col[j] = input[j][i];
      }

      this->insert({this->props[i], std::move(curr_col)});
    }
  } else {
    if (this->props.size() != input.size()) {
      throw std::runtime_error("Prop vector shall name all elements.");
    }

    auto name_it = this->props.begin();

    for (const auto &in_vec : input) {
      if (in_vec.size() != req_vec_size) {
        throw std::runtime_error(
            "Input vector doesn't match expected vector size.");
      }
      this->insert({*(name_it++), in_vec});
    }
  }
}

poet::Field::Field(std::uint32_t vec_s, const std::vector<double> &input,
                   const std::vector<std::string> &prop_vec)
    : Field(vec_s) {
  const uint32_t expected_size = prop_vec.size() * req_vec_size;

  if (expected_size != input.size()) {
    throw std::runtime_error(
        "Input vector have more (or less) elements than expected.");
  }

  auto name_it = prop_vec.begin();

  for (uint32_t i = 0; i < expected_size; i += req_vec_size) {
    auto input_pair = std::make_pair(
        *(name_it++), std::vector<double>(&input[i], &input[i + req_vec_size]));
    this->insert(input_pair);
  }

  this->props = prop_vec;
}

auto poet::Field::AsVector() const -> poet::FieldColumn {
  const uint32_t min_size = req_vec_size * this->size();

  poet::FieldColumn output;
  output.reserve(min_size);

  for (const auto &elem : props) {
    const auto map_it = this->find(elem);

    const auto start = map_it->second.begin();
    const auto end = map_it->second.end();

    output.insert(output.end(), start, end);
  }

  return output;
}

void poet::Field::update(const poet::Field &input) {
  for (const auto &input_elem : input) {
    auto it_self = this->find(input_elem.first);
    if (it_self == this->end()) {
      continue;
    }

    it_self->second = input_elem.second;
  }
}

auto poet::Field::As2DVector() const -> std::vector<poet::FieldColumn> {
  std::vector<poet::FieldColumn> output;
  output.reserve(this->size());

  for (const auto &elem : props) {
    const auto map_it = this->find(elem);
    output.push_back(map_it->second);
  }

  return output;
}

poet::FieldColumn &poet::Field::operator[](const std::string &key) {
  if (this->find(key) == this->end()) {
    props.push_back(key);
  }

  return std::unordered_map<std::string, FieldColumn>::operator[](key);
}

SEXP poet::Field::asSEXP() const {
  const std::size_t cols = this->props.size();

  SEXP s_names = PROTECT(Rf_allocVector(STRSXP, cols));
  SEXP s_output = PROTECT(Rf_allocVector(VECSXP, cols));

  for (std::size_t prop_i = 0; prop_i < this->props.size(); prop_i++) {
    const auto &name = this->props[prop_i];
    SEXP s_values = PROTECT(Rf_allocVector(REALSXP, this->req_vec_size));
    const auto values = this->find(name)->second;

    SET_STRING_ELT(s_names, prop_i, Rf_mkChar(name.c_str()));

    for (std::size_t i = 0; i < this->req_vec_size; i++) {
      REAL(s_values)[i] = values[i];
    }

    SET_VECTOR_ELT(s_output, prop_i, s_values);

    UNPROTECT(1);
  }

  SEXP s_rownames = PROTECT(Rf_allocVector(INTSXP, this->req_vec_size));
  for (std::size_t i = 0; i < this->req_vec_size; i++) {
    INTEGER(s_rownames)[i] = static_cast<int>(i + 1);
  }

  Rf_setAttrib(s_output, R_ClassSymbol,
               Rf_ScalarString(Rf_mkChar("data.frame")));
  Rf_setAttrib(s_output, R_NamesSymbol, s_names);
  Rf_setAttrib(s_output, R_RowNamesSymbol, s_rownames);

  UNPROTECT(3);

  return s_output;
}

poet::Field &poet::Field::operator=(const FieldColumn &cont_field) {
  if (cont_field.size() != this->size() * req_vec_size) {
    throw std::runtime_error(
        "Field::SetFromVector: vector does not match expected size");
  }

  uint32_t vec_p = 0;
  for (const auto &elem : props) {
    const auto start = cont_field.begin() + vec_p;
    const auto end = start + req_vec_size;

    const auto map_it = this->find(elem);

    map_it->second = FieldColumn(start, end);

    vec_p += req_vec_size;
  }

  return *this;
}

poet::Field &
poet::Field::operator=(const std::vector<FieldColumn> &cont_field) {
  if (cont_field.size() != this->size()) {
    throw std::runtime_error(
        "Input field contains more or less elements than this container.");
  }

  auto in_field_it = cont_field.begin();

  for (const auto &elem : props) {
    if (in_field_it->size() != req_vec_size) {
      throw std::runtime_error(
          "One vector contains more or less elements than expected.");
    }

    const auto map_it = this->find(elem);

    map_it->second = *(in_field_it++);
  }
  return *this;
}

void poet::Field::fromSEXP(const SEXP &s_rhs) {
  this->clear();

  SEXP s_vec = PROTECT(Rf_coerceVector(s_rhs, VECSXP));
  SEXP s_names = PROTECT(Rf_getAttrib(s_vec, R_NamesSymbol));

  std::size_t cols = static_cast<std::size_t>(Rf_length(s_vec));

  this->props.clear();
  this->props.reserve(cols);

  for (std::size_t i = 0; i < cols; i++) {
    const std::string prop_name(CHAR(STRING_ELT(s_names, i)));
    this->props.push_back(prop_name);

    SEXP s_values = PROTECT(VECTOR_ELT(s_vec, i));

    if (i == 0) {
      this->req_vec_size = static_cast<std::uint32_t>(Rf_length(s_values));
    }

    FieldColumn input(this->req_vec_size);

    for (std::size_t j = 0; j < this->req_vec_size; j++) {
      input[j] = static_cast<double>(REAL(s_values)[j]);
    }

    UNPROTECT(1);

    this->insert({prop_name, input});
  }

  UNPROTECT(2);
}
