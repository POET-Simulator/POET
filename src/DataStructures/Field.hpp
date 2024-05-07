#pragma once

#include <Rcpp.h>
#include <Rcpp/vector/Matrix.h>
#include <Rcpp/vector/instantiation.h>
#include <Rinternals.h>
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <iterator>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "DataStructures/PinnableMemory.hpp"

namespace poet {
template <class T> class Field {
public:
  Field() = delete;
  Field(Field &&) = delete;
  Field &operator=(Field &&) = delete;

  Field(const SEXP &s_rhs) {
    Rcpp::List in_list;
    try {
      in_list = Rcpp::List(s_rhs);
    } catch (const Rcpp::not_compatible &) {
      throw std::runtime_error("Input is not a list.");
    }

    if (in_list.size() == 0) {
      throw std::runtime_error("Input list is empty.");
    }

    Rcpp::NumericVector first_vec = Rcpp::NumericVector(in_list[0]);
    _m_vec_size = first_vec.size();

    // check if all elements are vectors of the same size
    for (const Rcpp::NumericVector &vec : in_list) {
      if (vec.size() != _m_vec_size) {
        throw std::runtime_error("Input vectors have different sizes.");
      }
    }

    _m_data.emplace(in_list.size() * _m_vec_size);

    // copy data
    for (std::size_t i = 0; i < in_list.size(); i++) {
      Rcpp::NumericVector curr_vec = Rcpp::NumericVector(in_list[i]);
      for (std::size_t j = 0; j < _m_vec_size; j++) {
        _m_data.value()[i * _m_vec_size + j] = curr_vec[j];
      }
    }

    // get the column names
    Rcpp::CharacterVector colnames = in_list.names();

    if (colnames.size() != in_list.size()) {
      throw std::runtime_error(
          "Column names do not match the number of columns.");
    }

    _m_colnames.resize(colnames.size());
    for (std::size_t i = 0; i < colnames.size(); i++) {
      _m_colnames[i] = (Rcpp::as<std::string>(colnames[i]));
    }
  }

  Field &operator=(Field &rhs) {
    if (this == &rhs) {
      return *this;
    }

    if (_m_vec_size != rhs._m_vec_size) {
      throw std::runtime_error("Vector sizes do not match.");
    }

    for (const auto &colname : _m_colnames) {
      auto it =
          std::find(rhs._m_colnames.begin(), rhs._m_colnames.end(), colname);

      if (it == rhs._m_colnames.end()) {
        continue;
      }

      std::span<T> rhs_col = rhs[colname];
      std::span<T> this_col = this->operator[](colname);

      std::memcpy(this_col.data(), rhs_col.data(), rhs_col.size_bytes());
    }

    return *this;
  }

  Field(Field &rhs)
      : _m_vec_size(rhs._m_vec_size), _m_colnames(rhs._m_colnames) {
    _m_data.emplace(rhs._m_data.value().size());
    std::memcpy(_m_data.value().data(), rhs._m_data.value().data(),
                rhs._m_data.value().size_bytes());
  }

  std::span<T> operator[](const std::string &colname) {
    auto it = std::find(_m_colnames.begin(), _m_colnames.end(), colname);
    if (it == _m_colnames.end()) {
      throw std::runtime_error("Column name not found.");
    }

    const std::size_t col_idx = std::distance(_m_colnames.begin(), it);

    return this->operator[](col_idx);
  }

  std::span<T> operator[](std::size_t col_idx) {
    if (col_idx >= _m_colnames.size()) {
      throw std::runtime_error("Column index out of bounds.");
    }

    return std::span<T>(&_m_data.value()[col_idx * _m_vec_size], _m_vec_size);
  }

  operator std::span<T>() { return std::span<T>(_m_data.value()); }

  operator SEXP() {
    Rcpp::List out_list(_m_colnames.size());

    for (std::size_t i = 0; i < _m_colnames.size(); i++) {
      Rcpp::NumericVector vec(_m_vec_size);
      std::memcpy(vec.begin(), _m_data.value().data() + i * _m_vec_size,
                  _m_vec_size * sizeof(T));
      out_list[i] = vec;
    }

    out_list.names() =
        Rcpp::CharacterVector(_m_colnames.begin(), _m_colnames.end());

    return std::move(out_list);
  }

  std::size_t size() const { return _m_colnames.size(); }

  std::size_t rows() const { return _m_vec_size; }

  T *data() { return _m_data.value().data(); }

  const std::vector<std::string> &colnames() const { return _m_colnames; }

private:
  std::optional<PinnableMemory<T>> _m_data;
  std::uint32_t _m_vec_size;
  std::vector<std::string> _m_colnames;
};
} // namespace poet