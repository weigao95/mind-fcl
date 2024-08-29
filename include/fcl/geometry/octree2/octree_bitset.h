//
// Created by Wei Gao on 2024/3/23.
//

#pragma once

#include <cstdint>

#if __cplusplus >= 202002L
#include <bit>
#endif

namespace fcl {
namespace octree2 {

/// Bitset specialized to a byte, used to represent the status of octree
/// children of a OctreeNode.
struct Bitset8 {
  static constexpr std::uint8_t bit_uint8 = 1;
  static constexpr std::uint8_t all_set = 0xFF;
  static constexpr std::uint8_t all_clear = 0x00;
  std::uint8_t bitset{all_clear};

  // Bit-wise method
  inline void set_i(std::uint8_t i_element) {
    bitset |= (bit_uint8 << i_element);
  }
  inline void clear_i(std::uint8_t i_element) {
    bitset &= (~(bit_uint8 << i_element));
  }
  inline bool test_i(std::uint8_t i_element) const {
    return static_cast<bool>(bitset & bit_uint8 << i_element);
  }

  // Count the #of bits in the element
  inline int count_of_bits() const {
#if __cplusplus >= 202002L
    return std::popcount(bitset);
#else
    // https://github.com/KAlO2/blog/blob/master/popcount/popcount.cpp
    const std::uint32_t n = bitset;
    uint32_t tmp = n - ((n >> 1) & 033333333333) - ((n >> 2) & 011111111111);
    return ((tmp + (tmp >> 3)) & 030707070707) % 63;
#endif
  }

  // Mutate and query as whole
  inline bool is_all_cleared() const { return bitset == all_clear; }
  inline bool is_all_set() const { return bitset == all_set; }
  inline void clear_all() { bitset = all_clear; }
  inline void set_all() { bitset = all_set; }

  // Constructor
  explicit Bitset8() = default;
  explicit Bitset8(std::uint8_t raw_bitset) : bitset(raw_bitset) {}
};

}  // namespace octree2
}  // namespace fcl
