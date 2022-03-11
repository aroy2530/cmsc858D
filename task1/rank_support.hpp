#ifndef RANK_SUPPORT_H
#define RANK_SUPPORT_H

#include "../include/compact_vector.hpp"
#include <iostream>
#include <string>

class rank_support {
    private: 
        uint64_t s; // size of superblock
        uint64_t b; // size of block
        uint64_t num_sb; // number of superblocks
        uint64_t num_b; // number of blocks
        compact::vector<uint64_t> r_s; // superblocks table
        compact::vector<uint64_t> r_b; // blocks table
        uint64_t read_int(uint64_t* word, uint8_t offset, uint8_t len);

    public:
        explicit rank_support(compact::vector<uint64_t,1>& new_bv); // constructor
        uint64_t rank1(uint64_t i);
        uint64_t overhead();
        uint64_t save(std::string& fname);
        uint64_t load(std::string& fname);
        compact::vector<uint64_t,1>& bv; // pointer to underlying bit-vector
        uint64_t n; // length of bit-vector

};

#endif