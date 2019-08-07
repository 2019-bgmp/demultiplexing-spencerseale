#!/usr/bin/env python3

# lines 1452986940
# reads 363246735

# Use this later for demultiplexing
#IDX_FILE = "/projects/bgmp/shared/2017_sequencing/indexes.txt"

#actual files
R1 = "/projects/bgmp/sseale/projects/demult/1294_S1_L008_R1_001.fastq.gz"
R2 = "/projects/bgmp/sseale/projects/demult/1294_S1_L008_R2_001.fastq.gz"
R3 = "/projects/bgmp/sseale/projects/demult/1294_S1_L008_R3_001.fastq.gz"
R4 = "/projects/bgmp/sseale/projects/demult/1294_S1_L008_R4_001.fastq.gz"

#test files
# R1 = "/projects/bgmp/sseale/projects/demult/r1.fq.gz"
# R2 = "/projects/bgmp/sseale/projects/demult/r2.fq.gz"
# R3 = "/projects/bgmp/sseale/projects/demult/r3.fq.gz"
# R4 = "/projects/bgmp/sseale/projects/demult/r4.fq.gz"

import numpy as np
import gzip
import matplotlib.pyplot as plt

# creating function to determine phred score from Q score
def convert_phred(letter):
    return ord(letter) - 33

# initializing an empty array to store values in

with gzip.open(R1, "rt") as r1, gzip.open(R2, "rt") as r2, gzip.open(R3, "rt") as r3, gzip.open(R4, "rt") as r4:
    LN1 = 0
    LN2 = 0
    LN3 = 0
    LN4 = 0
# creating array storing sum of phred scores for each pos in read 1 file
    for line in r1:
        # state where in the analysis it is at
        if LN1 % 1000000 == 0:
            print(R1 + ':' + str(LN1),flush=True)
        LN1 += 1
        if LN1 % 4 == 0:
            line = line.strip()
# When the qscore line is hit the first time, initialize the array, but only the first time, so LN = 4
            if LN1 == 4:
                array_len1 = len(line)
                phred_array_r1 = np.zeros((array_len1),dtype=float)
# Once the array is initialized, so for every qscore after the first, sum to array
            count = 0
            for x in line:
                phred_array_r1[count] += convert_phred(x)
                count += 1

# repeat for read 2, see above comments
    for line in r2:
        # state where in the analysis it is at
        if LN2 % 1000000 == 0:
            print(R2 + ':' + str(LN2),flush=True)
        LN2 += 1
        if LN2 % 4 == 0:
            line = line.strip()
            if LN2 == 4:
                array_len2 = len(line)
                phred_array_r2 = np.zeros((array_len2),dtype=float)
            count = 0
            for x in line:
                phred_array_r2[count] += convert_phred(x)
                count += 1

# repeat for read 3, see above comments
    for line in r3:
        # state where in the analysis it is at
        if LN3 % 1000000 == 0:
            print(R3 + ':' + str(LN3),flush=True)
        LN3 += 1
        if LN3 % 4 == 0:
            line = line.strip()
            if LN3 == 4:
                array_len3 = len(line)
                phred_array_r3 = np.zeros((array_len3),dtype=float)
            count = 0
            for x in line:
                phred_array_r3[count] += convert_phred(x)
                count += 1

# repeat for read 4, see above comments
    for line in r4:
        # state where in the analysis it is at
        if LN4 % 1000000 == 0:
            print(R4 + ':' + str(LN4),flush=True)
        LN4 += 1
        if LN4 % 4 == 0:
            line = line.strip()
            if LN4 == 4:
                array_len4 = len(line)
                phred_array_r4 = np.zeros((array_len4),dtype=float)
            count = 0
            for x in line:
                phred_array_r4[count] += convert_phred(x)
                count += 1

## testing to see whether arrays are made correctly having summed up values
# print(phred_array_r1)
# print(phred_array_r2)
# print(phred_array_r3)
# print(phred_array_r4)

# initializing arrays to store the mean in
mean_r1 = []
mean_r2 = []
mean_r3 = []
mean_r4 = []

# taking the average of the arrays created earlier that store total qscore for each read and taking average at each index
for _ in phred_array_r1:
    mean_r1.append(_/(LN1/4))
for _ in phred_array_r2:
    mean_r2.append(_/(LN2/4))
for _ in phred_array_r3:
    mean_r3.append(_/(LN3/4))
for _ in phred_array_r4:
    mean_r4.append(_/(LN4/4))

# print(mean_r1)
# print(mean_r2)
# print(mean_r3)
# print(mean_r4)

# plot read 1
plt.plot(list(range(array_len1)), mean_r1, "o")
plt.title("Read 1: Mean Q Score vs BP Number")
plt.xlabel("BP Number")
plt.ylabel("Mean Q Score")
plt.savefig("R1_mean_dist_plt.png")
plt.close()

# plot read 2
plt.plot(list(range(array_len2)), mean_r2, "o")
plt.title("Read 2: Mean Q Score vs BP Number")
plt.xlabel("BP Number")
plt.ylabel("Mean Q Score")
plt.savefig("R2_mean_dist_plt.png")
plt.close()

# plot read 3
plt.plot(list(range(array_len3)), mean_r3, "o")
plt.title("Read 3: Mean Q Score vs BP Number")
plt.xlabel("BP Number")
plt.ylabel("Mean Q Score")
plt.savefig("R3_mean_dist_plt.png")
plt.close()

# plot read 4
plt.plot(list(range(array_len4)), mean_r4, "o")
plt.title("Read 4: Mean Q Score vs BP Number")
plt.xlabel("BP Number")
plt.ylabel("Mean Q Score")
plt.savefig("R4_mean_dist_plt.png")
plt.close()

# # adding indices to dictionary, use this in fall when running demult with actual code
# idx_dict = {}
#
# with open(IDX_FILE, "r") as idxop:
#     LN = 0
#     for line in idxop:
#         line = line.strip()
#         LN += 1
#         if LN > 1:
#             parts = line.split()
#             idx_dict[parts[3]] = parts[4]
#
# #print(idx_dict)
