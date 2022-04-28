import sys
import pickle
import numpy as np
import fastaparser
import time

def main():
    index = sys.argv[1]
    queries_path = sys.argv[2]
    query_mode = sys.argv[3]
    output = sys.argv[4]

    prefs = []
    ranges = []

    # read in from index
    with open(index, 'rb') as f:
        ref = pickle.load(f)
        sa = np.load(f)
        k = pickle.load(f)
        if k != 0:
            prefs = np.load(f)
            ranges = np.load(f)

    # parser = ["aa", "ab", "aba", "b", "ba", "abaa"]
    # read in reference file
    with open(queries_path) as queries_file:
        parser = fastaparser.Reader(queries_file)

        num_queries = 1000
        query_counter = 0
        start_time = time.time()

        for seq in parser:
            if query_counter == num_queries:
                break
            id = seq.id
            query = seq.sequence_as_string()
            # id = "test"

            # list where output per query is stored
            q_output = [id]

            # query = seq
            # print(query)

            query_counter += 1

            # naive
            if query_mode == "naive":
                # set original l and r
                og_l = 0
                og_r = len(ref)-1

                # for testing
                # query = "ba"
                # k=0
                q_output = naive_search(ref, k, query, q_output, prefs, ranges, sa, og_l, og_r)

                # q_output = str(q_output)
                # print(q_output)

                file=open(output,'a')
                for x in q_output:
                    xs = str(x)
                    file.writelines(xs +'\t')
                file.write('\n')
                file.close()
                      
            
            elif query_mode == "simpaccel":
                print("not implemented")

    end_time = time.time()
    time_elapsed = end_time - start_time
    print("Time: ", time_elapsed)
            
                        
def naive_search(ref, k, query, q_output, prefs, ranges, sa, og_l, og_r):
    # prefix table check
    if k != 0:
        prefix = query[:k]
        # print(prefix)

        if len(query) < k:
            # not in prefix table and not in suffix array
            q_output.append(0)
            return q_output

        else:
            # find the index in the prefix table
            if prefix in prefs:
                preftab_idx = np.where(prefs==prefix)

                sa_range = ranges[preftab_idx][0]
                og_l = sa_range[0]
                og_r = sa_range[1]

            else:
                # a prefix of length k does not exist in the prefix table and thus the suffix array
                q_output.append(0)
                return q_output

    # print(og_l)
    # print(og_r)

    # vars  
    l = og_l
    r = og_r
    lower_bound = l # will mark first hit
    found = False # keep track of whether the query is present
    
    ## edge cases
    # if r - 1 is 0, only one element in the array
    if r - l == 0:
        if ref[sa[l]:].startswith(query):
            q_output.append(1)
            q_output.append(sa[l])
        else:
            q_output.append(0)
        return q_output

    # lower bound
    while l < r:
        c = int((l+r)/2)
        if query + "#" < ref[sa[c]:]:
            # c == og_l to avoid going over bounds
            # if the suffix at c-1 contains the query, we have not found the lower bound
            if c == og_l or (ref[sa[c]:].startswith(query) and not ref[sa[c-1]:].startswith(query)):
                lower_bound = c
                found = True
                break
            r = c
        else:
            if r - l == 1:
                if ref[sa[r]:].startswith(query):
                    lower_bound = r
                    found = True
                    break
                else:
                    break
            l = c

    if not found:
        q_output.append(0)
        return q_output

    # upper bound
    upper_bound = 0
    if found:
        # one occurrence
        if og_r - lower_bound == 0:
            q_output.append(1)
            q_output.append(sa[lower_bound])
            return q_output
        else:
            l = lower_bound
            r = og_r
            while l < r:
                c = int((l+r)/2)
                if query + "{" < ref[sa[c]:]:
                    r = c
                else:
                    if c == og_r or (ref[sa[c]:].startswith(query) and not ref[sa[c+1]:].startswith(query)):
                        upper_bound = c
                        break
                    if r - l == 1:
                        if ref[sa[r]:].startswith(query):
                            upper_bound = r
                            break
                        else:
                            break
                    l = c

            q_output.append(upper_bound - lower_bound + 1)
            
            # get list of hits
            for x in range(lower_bound, upper_bound + 1):
                q_output.append(sa[x])

    return q_output

# calculate least common prefix
# def LCP(a, b):
#     min_len = min(len(a), len(b))
#     for x in range(min_len):
#         if a[x] != b[x]:
#             return x
#     return min_len

# def simpaccel_binary_search(ref, sa, query, l, r):
#     # compute lcp between query and left bound
#     x = LCP(query, ref[sa[l]:])
#     y = LCP(query, ref[sa[r]:])

#     while l < r:
#         if query < ref[sa[c]:]:
#             if c == l + 1:
#                 return c
#             else:
#                 r = c
        
#         if query > ref[sa[c]:]:
#             if c == r - 1:
#                 return c
#             else:
#                 l = c

#         c = math.floor((l + r) / 2)

if __name__ == "__main__":
    main()