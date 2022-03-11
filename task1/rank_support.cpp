
#include "rank_support.hpp"
#include "../include/compact_vector.hpp"
#include <iostream>
#include <cmath>
#include <cstdint>

using namespace std;
// Implement a succinct, constant-time, bit-vector rank operation

namespace sdsl {
    const uint64_t lo_set[] = {0x0000000000000000ULL, 0x0000000000000001ULL, 0x0000000000000003ULL, 0x0000000000000007ULL, 0x000000000000000FULL, 0x000000000000001FULL, 0x000000000000003FULL, 0x000000000000007FULL, 0x00000000000000FFULL, 0x00000000000001FFULL, 0x00000000000003FFULL, 0x00000000000007FFULL, 0x0000000000000FFFULL, 0x0000000000001FFFULL, 0x0000000000003FFFULL, 0x0000000000007FFFULL, 0x000000000000FFFFULL, 0x000000000001FFFFULL, 0x000000000003FFFFULL, 0x000000000007FFFFULL, 0x00000000000FFFFFULL, 0x00000000001FFFFFULL, 0x00000000003FFFFFULL, 0x00000000007FFFFFULL, 0x0000000000FFFFFFULL, 0x0000000001FFFFFFULL, 0x0000000003FFFFFFULL, 0x0000000007FFFFFFULL, 0x000000000FFFFFFFULL, 0x000000001FFFFFFFULL, 0x000000003FFFFFFFULL, 0x000000007FFFFFFFULL, 0x00000000FFFFFFFFULL, 0x00000001FFFFFFFFULL, 0x00000003FFFFFFFFULL, 0x00000007FFFFFFFFULL, 0x0000000FFFFFFFFFULL, 0x0000001FFFFFFFFFULL, 0x0000003FFFFFFFFFULL, 0x0000007FFFFFFFFFULL, 0x000000FFFFFFFFFFULL, 0x000001FFFFFFFFFFULL, 0x000003FFFFFFFFFFULL, 0x000007FFFFFFFFFFULL, 0x00000FFFFFFFFFFFULL, 0x00001FFFFFFFFFFFULL, 0x00003FFFFFFFFFFFULL, 0x00007FFFFFFFFFFFULL, 0x0000FFFFFFFFFFFFULL, 0x0001FFFFFFFFFFFFULL, 0x0003FFFFFFFFFFFFULL, 0x0007FFFFFFFFFFFFULL, 0x000FFFFFFFFFFFFFULL, 0x001FFFFFFFFFFFFFULL, 0x003FFFFFFFFFFFFFULL, 0x007FFFFFFFFFFFFFULL, 0x00FFFFFFFFFFFFFFULL, 0x01FFFFFFFFFFFFFFULL, 0x03FFFFFFFFFFFFFFULL, 0x07FFFFFFFFFFFFFFULL, 0x0FFFFFFFFFFFFFFFULL, 0x1FFFFFFFFFFFFFFFULL, 0x3FFFFFFFFFFFFFFFULL, 0x7FFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL};
}

// constructor
rank_support::rank_support(compact::vector<uint64_t,1>& new_bv) 
    : bv(new_bv)
    , n(bv.size())
    , s(pow(ceil(log2(n)), 2)) // s = ceiling(log 2 (n))
    , num_sb(ceil(static_cast< float >(n)/static_cast< float >(s)))
    , r_s(ceil(log2(n)), num_sb) // r_s counter width = ceil(log(n))
    , b(ceil(log2(n)))
    , num_b(ceil(static_cast< float >(n)/static_cast< float >(b)))
    , r_b(ceil(log2(pow(log2(n), 2))), num_b)  // r_b cw = ceil(log(log2(n)))
 {

    // print out size of r_s to make sure it worked
    // cout << "\nr_s size: " << r_s.size(); //1 1
    // cout << "\nr_s bits: " << r_s.bits(); //4 5
    // cout << "\nr_b size: " << r_b.size(); //4 4
    // cout << "\nr_b bits: " << r_b.bits(); //4 5

    // cout << "\nprinting r_s\n";
    // for (int i = 0; i < r_s.size(); i++) {
    //     cout << r_s.at(i) << " ";
    // }

    uint64_t rank = 0;  // keep track of ongoing rank to avoid recalculating
    int curr_sb_idx = 0; // keep track of starting index of current superblock

    // loop through length of bitvector
    for (int j = 0; j < n; j++) {

        // if j is the starting index of a superblock, assign rank
        if (j % s == 0) {
            r_s.at(j/s) = rank;
            curr_sb_idx = j/s;
        }

        // if j is the starting index of a block, assign (rank - rank at start of current superblock)
        if (j % b == 0) {
            r_b.at(j/b) = rank - r_s.at(curr_sb_idx);
        }

        // increment rank as necessary (done at the end to avoid messing up superblocks/blocks)
        if (bv.at(j) == 1) {
            rank++;
        }
    }

    // cout << "\nprinting r_s\n";
    // for (int i = 0; i < r_s.size(); i++) {
    //     cout << r_s.at(i) << " ";
    // }

    // cout << "\nprinting r_b\n";
    // for (int i = 0; i < r_b.size(); i++) {
    //     cout << r_b.at(i) << " ";
    // }

}

// Returns the number of 1s in the underlying bit-vector up to position i (inclusive)
uint64_t rank_support::rank1(uint64_t i) {
    uint64_t b_idx = i/b; // get starting index of block i falls in
    // cout << "\nb" << b_idx;
    uint64_t s_idx = i/s; // get starting index of superblock i falls in
    // cout << s_idx;

    uint64_t bv_offset = b_idx * b; // index of b_idx in respect to bv
    uint64_t bl = i - bv_offset + 1; // how much to slice 
    uint64_t* bv_data_ptr = bv.get();

    // taken from read_int function in sdsl-lite
    uint64_t res = read_int(bv_data_ptr+(bv_offset>>6), bv_offset&0x3F, bl);
    // cout << "\n" <<  __builtin_popcount(res);

    return r_s.at(s_idx) + r_b.at(b_idx) + __builtin_popcount(res);
}

uint64_t rank_support::overhead() {
    return (ceil(log2(n)) * num_sb) + (ceil(log2(pow(log2(n), 2))) * num_b) + n;
}

uint64_t rank_support::save(string& fname) {
    ofstream save_file(fname);

    bv.serialize(save_file);
    r_s.serialize(save_file);
    r_b.serialize(save_file);

    save_file.close();
}

uint64_t rank_support::load(string& fname) {
    bv.deserialize(fname);
    n = bv.size();

    r_s.deserialize(fname);
    s = pow(ceil(log2(n)), 2);
    num_sb = ceil(static_cast< float >(n)/static_cast< float >(s));

    r_b.deserialize(fname);
    b = ceil(log2(n));
    num_b = ceil(static_cast< float >(n)/static_cast< float >(b));
}

uint64_t rank_support::read_int(uint64_t* word, uint8_t offset, uint8_t len) {
    uint64_t w1 = (*word)>>offset;
    if ((offset+len) > 64) { // if offset+len > 64
        return w1 |  // w1 or w2 adepted:
               ((*(word+1) & sdsl::lo_set[(offset+len)&0x3F])   // set higher bits zero
                << (64-offset));  // move bits to the left
    } else {
        return w1 & sdsl::lo_set[len];
    }
}

int main() {
    // compact::vector<uint64_t, 1> bb(16); // same as
    // compact::vector<uint64_t, 1> bb();
    // // bb.resize(16);                    // same as
    // // compact::vector<uint64_t> bb(1, 16);

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

    for (int i = 0; i < b.size(); i++) {
        cout << b[i];
    }

    rank_support r(b); // let r access b via a pointer
    auto x = r.rank1(8); // x now holds the value of rank1(b, 1);
                   
    cout << "\ncount " << x;
}