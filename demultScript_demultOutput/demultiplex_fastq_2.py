#!/usr/bin/env python3

import os
import gzip
import numpy as np
import argparse

#function to take in input when script is called
def get_args():
    parser = argparse.ArgumentParser(description="demultiplex and quality filter fastq files")
    parser.add_argument("-d", "--dir", help="Input directory containing all fastq files to be read", type=str)
    parser.add_argument("-r1", "--read1", help="Input r1 fastq file", type=str)
    parser.add_argument("-r2", "--read2", help="Input r2 fastq file", type=str)
    parser.add_argument("-r3", "--read3", help="Input r3 fastq file", type=str)
    parser.add_argument("-r4", "--read4", help="Input r4 fastq file", type=str)
    return parser.parse_args()
args = get_args()

#setting coverage cutoff
#if mean of r1/r4 is at or above the bio_cutoff, quality passes
#if single nt in r2/r3 is below idx_cutoff, quality fails
idx_cutoff = 30
bio_cutoff = 25

#file containing indices used in sequencing run
IDX_FILE = "/projects/bgmp/shared/2017_sequencing/indexes.txt"

#actual files
# R1 = "1294_S1_L008_R1_001.fastq.gz"
# R2 = "1294_S1_L008_R2_001.fastq.gz"
# R3 = "1294_S1_L008_R3_001.fastq.gz"
# R4 = "1294_S1_L008_R4_001.fastq.gz"
# direct = "/projects/bgmp/sseale/projects/demult/"

R1 = args.read1
R2 = args.read2
R3 = args.read3
R4 = args.read4
direct = args.dir

#test files
# R1 = "/projects/bgmp/sseale/projects/demult/demultiplexing-spencerseale/unit-tests/r1.fq.gz"
# R2 = "/projects/bgmp/sseale/projects/demult/demultiplexing-spencerseale/unit-tests/r2.fq.gz"
# R3 = "/projects/bgmp/sseale/projects/demult/demultiplexing-spencerseale/unit-tests/r3.fq.gz"
# R4 = "/projects/bgmp/sseale/projects/demult/demultiplexing-spencerseale/unit-tests/r4.fq.gz"

#dictionary for nucleotide compliments
seqKey = {
"A": "T",
"T": "A",
"G": "C",
"C": "G",
"N": "N"
}

#function to read record
def read_record(read_file, read_list):
    read_list.append(read_file.readline().strip())
    read_list.append(read_file.readline().strip())
    read_list.append(read_file.readline().strip())
    read_list.append(read_file.readline().strip())
    return read_list

#function to write out corrected read to file
def write_record(file_r1, file_r2):
    file_r1.write(cor_r1list[0])
    file_r1.write("\n")
    file_r1.write(cor_r1list[1])
    file_r1.write("\n")
    file_r1.write(cor_r1list[2])
    file_r1.write("\n")
    file_r1.write(cor_r1list[3])
    file_r1.write("\n")
    file_r2.write(cor_r4list[0])
    file_r2.write("\n")
    file_r2.write(cor_r4list[1])
    file_r2.write("\n")
    file_r2.write(cor_r4list[2])
    file_r2.write("\n")
    file_r2.write(cor_r4list[3])
    file_r2.write("\n")

#function will take in a sequence and return its reverse compliment
def rev_comp(seq):
    revComp = ""
    for nt in range(len(seq)-1, -1, -1):
        revComp += seqKey[seq[nt]]
    return revComp

#function to determine phred score from Q score
def convert_phred(letter):
    return ord(letter) - 33

#function will take the two biological read records currently being looked at and stored
#in memory and append the index 1 sequence as appears and the rev comp of index 2 seq to
#header of each biological read separated by a "-" (i.e. AAAAAAAA-AAAAAAAA)
def str_add(r1memory, r4memory):
    cor_r1list = []
    cor_r4list = []
    idx_format = "_"+r2list[1]+"-"+rev_comp(r3list[1])
    cor_r1list.append(r1list[0]+idx_format)
    cor_r1list.append(r1list[1])
    cor_r1list.append(r1list[2])
    cor_r1list.append(r1list[3])
    cor_r4list.append(r4list[0]+idx_format)
    cor_r4list.append(r4list[1])
    cor_r4list.append(r4list[2])
    cor_r4list.append(r4list[3])
    return cor_r1list, cor_r4list

#function to build an array to hold quality scores
def array_builder(list):
    array_length_variable = len(list[3])
    array = np.zeros((array_length_variable),dtype=float)
    return array

#grabbing indexes used in seq from separate file and importing into dictionary to allow for fast index look up compared to list
#keys are index seq and values are index seq
idx_dict = {}
with open(IDX_FILE, "r") as idxop:
   LN = 0
   for line in idxop:
       line = line.strip()
       LN += 1
       if LN > 1:
           parts = line.split()
           idx_dict[parts[4]] = parts[4]

#creating all the files from the indices as well as the unkown and index hop files
os.chdir("/projects/bgmp/sseale/projects/demult/demult_files")
unknown_r1 = open("unknown_low-quality_r1.fq", "w")
unknown_r2 = open("unknown_low-quality_r2.fq", "w")
hop_r1 = open("idx_hop_r1.fq", "w")
hop_r2 = open("idx_hop_r2.fq", "w")

#this is where the index files are created. also adding file names to a dictionary
file_pointer = {}
for x in idx_dict.values():
    r1 = x+"_r1.fastq"
    r2 = x+"_r2.fastq"
    z = open(r1, "w")
    q = open(r2, "w")
    file_pointer[x] = (r1, r2)
    z.close()
    q.close()

#setting summary stats
num_reads_low_qual_unkwn = 0
num_reads_idx_hop = 0
num_reads_pass = 0

#initializing lists to store each read currently being analyzed to be stored in memory
r1list = []
r2list = []
r3list = []
r4list = []

#opening 4 specified input files for eading and running demultiplexing loop
with gzip.open(direct+R1, "rt") as r1op, gzip.open(direct+R2, "rt") as r2op, gzip.open(direct+R3, "rt") as r3op, gzip.open(direct+R4, "rt") as r4op:
    LN1 = 0
    LN2 = 0
    LN3 = 0
    LN4 = 0
    #loop will break when end of one of the files is reached
    while True:
        LN1 += 1
        LN2 += 1
        LN3 += 1
        LN4 += 1
        #reading r1
        r1list = read_record(r1op, r1list)
        if r1list[0] == "":
            break
        #must rebuild array each loop to handle different read lengths
        phred_array_r1 = array_builder(r1list)
        count = 0
        for qscore in r1list[3]:
            phred_array_r1[count] = convert_phred(qscore)
            count += 1
        #q score check, checking mean of all nt
        if np.mean(phred_array_r1) >= bio_cutoff:
            r1_low_qual = False
        else:
            r1_low_qual = True
        #reading r2
        r2list = read_record(r2op, r2list)
        phred_array_r2 = array_builder(r2list)
        count = 0
        for qscore in r2list[3]:
            phred_array_r2[count] = convert_phred(qscore)
        #q score check, checking each nt of index for passing q score
            if phred_array_r2[count] >= idx_cutoff:
                r2_low_qual = False
                count += 1
            else:
                r2_low_qual = True
                break
        #reading r3
        r3list = read_record(r3op, r3list)
        phred_array_r3 = array_builder(r3list)
        count = 0
        for qscore in r3list[3]:
            phred_array_r3[count] = convert_phred(qscore)
        #q score check, checking each nt of index for passing q score
            if phred_array_r3[count] >= idx_cutoff:
                r3_low_qual = False
                count += 1
            else:
                r3_low_qual = True
                break
        #reading r4
        r4list = read_record(r4op, r4list)
        phred_array_r4 = array_builder(r4list)
        count = 0
        for qscore in r4list[3]:
            phred_array_r4[count] = convert_phred(qscore)
            count += 1
        #q score check, checking mean of all nt
        if np.mean(phred_array_r4) >= bio_cutoff:
            r4_low_qual = False
        else:
            r4_low_qual = True

        #running str_add function to turn r1 and r4 lists into the corrected format
        cor_r1list, cor_r4list = str_add(r1list, r4list)
        #checking whether any of the reads are low quality, if any true, all reads written to low quality files
        if r1_low_qual | r2_low_qual | r3_low_qual | r4_low_qual:
            write_record(unknown_r1, unknown_r2)
            num_reads_low_qual_unkwn += 1
        #checking if each index is unknown
        elif rev_comp(r3list[1]) not in idx_dict:
            write_record(unknown_r1, unknown_r2)
            num_reads_low_qual_unkwn += 1
        elif r2list[1] not in idx_dict:
            write_record(unknown_r1, unknown_r2)
            num_reads_low_qual_unkwn += 1
        #now both r2 and r3 indices have been checked to be used in seq experiment
        #if this if statment passes, will write to appropriate index file
        elif r2list[1] == rev_comp(r3list[1]):
            #appending to the end of each file
            read1write = open(file_pointer[r2list[1]][0], "a")
            read2write = open(file_pointer[r2list[1]][1], "a")
            write_record(read1write, read2write)
            read1write.close()
            read2write.close()
            num_reads_pass += 1
        #if this statement passes, then both indices are valid, but index hopping has occured
        elif r2list[1] != rev_comp(r3list[1]):
            write_record(hop_r1, hop_r2)
            num_reads_idx_hop += 1
        #wiping whats in memory to start over for next record
        r1list = []
        r2list = []
        r3list = []
        r4list = []

unknown_r1.close()
unknown_r2.close()
hop_r1.close()
hop_r2.close()

#reporting processing summary
num_reads = LN1 - 1
print(f"Number of total read pairs processed: {num_reads}")
print(f"Number of read pairs passing quality check and written out: {num_reads_pass}, {num_reads_pass/num_reads*100}%")
print(f"Number of read pairs omitted due to index hopping: {num_reads_idx_hop}, {num_reads_idx_hop/num_reads*100}%")
print(f"Number of read pairs omitted due to unknown index seq or low quality: {num_reads_low_qual_unkwn}, {num_reads_low_qual_unkwn/num_reads*100}%\n")

#reporting percentages for each index written out
for pointer in file_pointer:
    line_count = 0
    with open(file_pointer[pointer][0]) as openfile:
        for line in openfile:
            line_count += 1
        percent = ((line_count/4)/num_reads_pass)*100
        print(f"Sample having index {pointer}, yielded {percent}% of the total read pairs passing QC and written out.")

print("\nAnalysis complete.")
