#ifndef DATASTRUCTURES_H_
#define DATASTRUCTURES_H_

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace poet {

template <typename value_type> class NamedVector {
public:
  void insert(std::pair<std::string, value_type> to_insert) {
    this->names.push_back(to_insert.first);
    this->values.push_back(to_insert.second);
    container_size++;
  }

  bool empty() const { return this->container_size == 0; }

  std::size_t size() const { return this->container_size; }

  value_type &operator[](std::size_t i) { return values[i]; }

  const std::vector<std::string> &getNames() const { return this->names; }
  const std::vector<value_type> &getValues() const { return this->values; }

private:
  std::size_t container_size{0};

  std::vector<std::string> names{};
  std::vector<value_type> values{};
};
} // namespace poet
#endif // DATASTRUCTURES_H_
