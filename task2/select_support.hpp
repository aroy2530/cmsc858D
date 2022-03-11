#ifndef SELECT_SUPPORT_H
#define SELECT_SUPPORT_H

#include "../include/compact_vector.hpp"
#include "../task1/rank_support.hpp"
#include <iostream>
#include <string>

class select_support {
    private: 
        rank_support& r; // pointer to rank

    public:
        select_support(rank_support& new_r); // constructor
        uint64_t select1(uint64_t i);
        uint64_t overhead();
        uint64_t save(std::string& fname);
        uint64_t load(std::string& fname);

};

#endif