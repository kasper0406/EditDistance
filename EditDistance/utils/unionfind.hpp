#pragma once

#include <memory>
#include <sstream>
#include <string>

#include "math.hpp"

namespace Utils {
  using namespace std;
  
  template <class T>
  class UnionFind
  {
  public:
    class Element {
    public:
      Element(T e) : element(e), parent(this), rank(0)
      { }
      
      T element;
      
    private:
      Element* parent;
      uint64_t rank;
      
      friend class UnionFind;
    };
    
    static Element* Find(Element* x) {
      if (x->parent != x)
        x->parent = Find(x->parent);
      return x->parent;
    }
    
    static Element* Union(Element* A, Element* B) {
      auto A_root = Find(A);
      auto B_root = Find(B);
      if (A_root == B_root)
        return A_root;
      
      if ((A_root->rank) < (B_root->rank)) {
        A_root->parent = B_root;
        return B_root;
      } else if ((A_root->rank) > (B_root->rank)) {
        B_root->parent = A_root;
        return A_root;
      } else {
        B_root->parent = A_root;
        A_root->rank++;
        return A_root;
      }
    }
  };
  
  class IntervalUnionFind {
  public:
    IntervalUnionFind(int64_t range) : num_blocks(Math::ceil_div<int64_t>(range, block_size))
    {
      // Initialize blocks
      blocks_.reserve(num_blocks);
      for (int64_t i = 0; i < num_blocks; ++i) {
        blocks_.push_back(Block(i));
      }
    }
    
    void Union(int64_t A, int64_t B) {
      assert(A < B);
      auto A_max = Find(A);
      auto B_max = Find(B);
      assert(A_max < B_max);
      
      uint64_t A_block, A_bit_index, B_block, B_bit_index;
      tie(A_block, A_bit_index, ignore) = decompose(A_max);
      tie(B_block, B_bit_index, ignore) = decompose(B_max);
      
      // Unset the bit
      blocks_[A_block].bits &= ~((uint64_t)1 << A_bit_index);
      if (A_block != B_block) {
        auto new_root = UF::Union(blocks_[A_block].right.get(), blocks_[B_block].left.get());
        new_root->element = max((uint64_t)(B_max / block_size), B_block);
      }
    }
    
    int64_t Find(int64_t x) {
      uint64_t block, bit_index, element_index;
      tie(block, bit_index, element_index) = decompose(x);
      
      // cout << print_bits(blocks_[block].bits << element_index) << endl;
      int nonzero = Math::msb(blocks_[block].bits << element_index);
      if (nonzero == -1) {
        // All entries in the block are zero
        // cout << "Quering to block: " << UF::Find(blocks_[block].right.get())->element << endl;
        return Find(UF::Find(blocks_[block].right.get())->element * block_size);
      }
      
      return (64 - nonzero) + element_index + block_size * block;
    }
    
    void print() {
      for (int64_t i = 0; i < num_blocks; ++i) {
        cout << print_bits(blocks_[i].bits);
        if (i != num_blocks - 1)
          cout << " | ";
      }
      cout << endl;
    }
    
  private:
    const int64_t num_blocks;
    
    tuple<int64_t, int64_t, int64_t> decompose(int64_t index) const {
      uint64_t block = index / block_size;
      uint64_t bit_index = 64 - (index % block_size) - 1;
      uint64_t element_index = index % block_size;
      return { block, bit_index, element_index };
    }
    
    static constexpr int8_t block_size = 64;
    
    typedef UnionFind<uint64_t> UF;
    struct Block {
      Block(int64_t block) : bits(~(((uint64_t)1 << (64 - block_size)) - (uint64_t)1)),
                             left(unique_ptr<UF::Element>(new UF::Element(block))),
                             right(unique_ptr<UF::Element>(new UF::Element(block)))
      { }
      
      uint64_t bits;
      unique_ptr<UF::Element> left;
      unique_ptr<UF::Element> right;
    };
    vector<Block> blocks_;
    
    string print_bits(uint64_t bits) {
      stringstream ss;
      for (int64_t i = 63; i > 63 - block_size; --i) {
        if ((bits >> i) % 2 == 0)
          ss << "0";
        else
          ss << "1";
      }
      return ss.str();
    }
  };
}
