# read in a genome
# build the suffix array
# write string and suffix array to binary file

import sys
import numpy as np
from pydivsufsort import divsufsort
import fastaparser
import time
import pickle
import os

preftab_flag = False # flag to keep track of whether building a prefix table or not
ref_idx = 1 # arg index of reference
k = 0

if (sys.argv[1] == "--preftab"):
    preftab_flag = True
    k = int(sys.argv[2])
    ref_idx = 3

ref_path = sys.argv[ref_idx]
output_file = sys.argv[ref_idx + 1]

###############

# read in reference file
with open(ref_path) as ref_file:
    parser = fastaparser.Reader(ref_file)
    for seq in parser:
        ref = seq.sequence_as_string() + "$" # append $ sentinel

print("read in ref")
# ref = ref[:2319838]
print("Length of reference: ", len(ref))
# ref = "abaaba$"

start_time = time.time()

# build suffix array
sa = divsufsort(ref)
print("built suffix array")
# print(sa)
# print(len(sa))

if preftab_flag:
    st_idx = 0 # start index of prefix range
    end_idx = 0 # end index of prefix range
    curr_pref = ""
    prefs = []
    ranges = []

    # iterate through indices of suffix array
    for sa_idx in range(len(sa)):
        # suffix = ref[sa[sa_idx]:] 
        # print(suffix)

        # only if length of ref - length of suffix is greater than k, then add
        if len(ref) - sa[sa_idx] > k:
            prefix = ref[sa[sa_idx]:(sa[sa_idx] + k)]  # calculate prefix
            # print(prefix)

            # if prefix is not the same as previous prefix
            if prefix != curr_pref:
                # store in prefix table
                if curr_pref != "":
                    ranges.append([st_idx, end_idx])
                    prefs.append(curr_pref)
                if sa_idx == sa.size - 1:
                    ranges.append([sa_idx, sa_idx])
                    prefs.append(prefix)
                curr_pref = prefix
                st_idx = sa_idx
                end_idx = sa_idx
            elif prefix == curr_pref:
                # print("same")
                end_idx += 1
                if sa_idx == sa.size - 1:
                    ranges.append([st_idx, end_idx])
                    prefs.append(prefix)
    
    # print(preftab)

    prefs_np = np.asarray(prefs)
    ranges_np = np.asarray(ranges)
    print("finished prefix table")
    # print(prefs_np)
    # print(ranges_np)
    # print(preftab_np.size)


# timing
end_time = time.time()
time_elapsed = end_time - start_time
print("Time: ", time_elapsed)

# write binary to output file
with open(output_file, "wb") as output:
    pickle.dump(ref, output)
    np.save(output, sa)
    pickle.dump(k, output)
    if preftab_flag:
        np.save(output, prefs)
        np.save(output, ranges)

# file size
print("File size: ", os.path.getsize(output_file)/1000000, " MB")
