#pragma once

#include <cstddef>
#include <memory>
#include <span>

/**
 * @brief Class to allocate memory which can be pinned (e.g. for RDMA access)
 *
 * @tparam T Type of the memory to allocate
 */
template <class T> class PinnableMemory {
public:
  PinnableMemory() = delete;
  PinnableMemory(PinnableMemory &&) = delete;
  PinnableMemory &operator=(const PinnableMemory &) = delete;
  PinnableMemory &operator=(PinnableMemory &&) = delete;
  PinnableMemory(const PinnableMemory &other) = delete;

  /**
   * @brief Construct a new pinnable memory object
   *
   * @param count Count of elements to allocate
   */
  PinnableMemory(std::size_t count) : _m_data(new T[count]), _m_count(count){};

  ~PinnableMemory() = default;

  /**
   * @brief Return the base address of the allocated memory
   *
   * @return void* Base address of the allocated memory.
   */
  T *data() { return _m_data.get(); }

  /**
   * @brief Returns the element at the given index
   *
   * @param index Index of the element to return
   * @return T& Element at the given index
   */
  T &operator[](std::size_t index) { return _m_data.get()[index]; }

  /**
   * @brief Return a span over allocated memory
   *
   * @return std::span<T> Span over allocated memory
   */
  operator std::span<T>() { return std::span<T>(_m_data.get(), _m_count); }

private:
  // byte addressable unique pointer
  std::unique_ptr<T[]> _m_data;

  std::size_t _m_count;
};