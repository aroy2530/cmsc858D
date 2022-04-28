#include "../include/compact_vector.hpp"
#include "../task1/rank_support.hpp"
#include "select_support.hpp"
#include <iostream>
#include <string>
#include <cmath>

// constructor
select_support::select_support(rank_support& new_r) : r(new_r) {
}

// binary search-esque method using rank
uint64_t select_support::select1(uint64_t i) {
    uint64_t begin = 0;
    uint64_t end = r.n;
    uint64_t middle = (r.n)/2;

    //  If rank(n/2) > i, then the i th 1-bit is in B[1..n/2]
    uint64_t curr_rank = r.rank1(middle);
    while (curr_rank != i && (begin != end)) {
        if (curr_rank > i) {
            end = middle - 1;   // look in upper half
        }
        else if (curr_rank < i) {
            begin = middle + 1; // look in lower half
        }

        middle = floor((begin + end) / 2);  // reset middle
    }
    return middle;
}

uint64_t select_support::overhead() {
    return r.overhead();
}

uint64_t select_support::save(std::string& fname) {
    return r.save(fname);
}

uint64_t select_support::load(std::string& fname) {
    return r.load(fname);
}

int main () {
    compact::vector<uint64_t, 1> b(16);

    // cout << "bit vector size " << b.size() << "\n";
    // cout << "bit vector bits " << b.bits() << "\n";
    
    for (int i = 0; i < b.size(); i++) {
        b[i] = 0;
    }

    b[0] = 1;
    b[3] =1;
    b[5]=1;
    b[6]=1;
    b[7]=1;
    b[9]=1;
    b[12]=1;
    b[14]=1;

    // for (int i = 0; i < b.size(); i++) {
    //     std::cout << b[i];
    // }

    rank_support r(b); // let r access b via a pointer
    select_support s(r);
    std::cout << s.select1(3);
}
