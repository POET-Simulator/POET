#include "DataStructures/PinnableMemory.hpp"

#include <doctest/doctest.h>

TEST_CASE("Test pinnable memory") {
  // Test case 1: Create PinnableMemory object with count = 5
  PinnableMemory<int> memory(5);

  // Test case 2: Verify the base address of the allocated memory
  CHECK_NE(memory.data(), nullptr);

  // Test case 3: Access and modify elements using the subscript operator
  memory[0] = 10;
  memory[1] = 20;
  memory[2] = 30;
  memory[3] = 40;
  memory[4] = 50;

  CHECK_EQ(memory[0], 10);
  CHECK_EQ(memory[1], 20);
  CHECK_EQ(memory[2], 30);
  CHECK_EQ(memory[3], 40);
  CHECK_EQ(memory[4], 50);

  // Test case 4: Convert PinnableMemory to std::span and iterate over elements
  std::span<int> span = memory;
  int sum = 0;
  for (const auto &element : span) {
    sum += element;
  }
  CHECK_EQ(sum, 150);
}
