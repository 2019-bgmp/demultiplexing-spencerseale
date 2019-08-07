**Part 1**

*1*

| File | Read |
| --- | --- |
| 1294_S1_L008_R1_001.fastq.gz | read1 |
| 1294_S1_L008_R2_001.fastq.gz | index1 |
| 1294_S1_L008_R3_001.fastq.gz | index2 |
| 1294_S1_L008_R4_001.fastq.gz | read2 |

*2b*

Proposed index read qscore cutoff: 30\n
-Picking a higher qscore cutoff for the index read since it is more critical to ensure that the index sequences are correct. Having an incorrect index sequence is important for multiplexing. As long as not many indexes have errors, it will be simple to distinguish the correct index from error indexes in the final collection of data. Additionally, there are less nucleotides in the index sequence and therefore a higher chance that sequencing errors will occur in the same basepair, resulting in more repeats of the same error.

Proposed biological read qscore cutoff: 25\n
-A lower quality score for the biological reads are okay because the sequence is of a longer length and therefore the chances of the same sequencing errors occurring at the same bp position are less likely. It will be less difficult to determine the reads that have sequencing errors in alignment of the biological reads.

*2c*

```
command: zcat /projects/bgmp/sseale/projects/demult/1294_S1_L008_R2_001.fastq.gz /projects/bgmp/sseale/projects/demult/1294_S1_L008_R3_001.fastq.gz | grep -B 1 "^+" | grep -v "^-" | grep -v "^+" | grep -c "N"

output: 7304664
```

**Part 2**

Problem defined concisely:\n
-Demultiplex the PE sequencing run, separating the perfect match read pairs, unknown/low quality read pairs, and index hopping read pairs. Calculate number of occurrences of each.

Problem defined in detail:\n
-For each Paired-end (PE) read, extract the 8 bp p7 and p5 index reads and append them to the header of both the biological read 1 and 2 in the format: "TAAAAAAA-TAAAAAAA"\n
-Check whether the index sequence of both PE reads corresponds to the sequences of the indexes used in library prep.\n
-Separate each dual matched index read into a file for both biological read 1 and read 2 for each sequence resulting in 48 combined read 1 and read 2 fastq files.\n
-If the indexes on biological read 1 and read 2 are not perfect matches, but both indexes are listed in the indexes used in the sequencing run, index hopping occurred. For each instance of index hopping, add the biological read 1 fastq record to a dedicated file for read 1 index hopping and the biological read 2 fastq record to a second dedicated file fore read 2 index hopping.\n
-If one or both indexes do not match the index sequences used in library prep, add the read 1 and read 2 fastq records to separate fastq files for unknown index sequence. Additionally, if one or both of the indexes do not meet the quality cutoff add both read 1 and read 2 to the same unknown index sequence fastq files. There should be one read 1 and one read 2 file.\n
-For the properly matched (per index combination), index hopping, and unknown/low quality fastq files (52 files total) calculate the number of sequencing read pairs fitting into each category.\n

Helpful output:\n
-Number of read pairs in each of the above pairs after parsing.\n
-Reverse compliment of read 2 index to match to read 1 index to determine perfect index match.\n
-Number of lines in each output file

*pseudocode*

-Initialize an empty list (r1list) that will hold the current r1 file record in memory.\n
-Initialize an empty list (r4list) that will hold the current r4 file record in memory.\n
-Initialize an empty list (r2list) that will hold the current r2 file record in memory.\n
-Initialize an empty list (r3list) that will hold the current r3 file record in memory.\n

-Create a dictionary for all indexes used in exp, with the index name as key, sequence as value. Create dictionary directly from provided index file.
```
idx_dict = {}

with open(IDX_FILE, "r") as idxop:
   LN = 0
   for line in idxop:
       line = line.strip()
       LN += 1
       if LN > 1:
           parts = line.split()
           idx_dict[parts[3]] = parts[4]
```

-Create the following files:\n
  -Using dictionary made above for the index sequences, create a r1 and r2 file for each index listed in given index file in the format (48 total): idx1seq_r1.fq and idx1seq_r2.fq.\n
    -idx1seq = actual sequence of idx 1, don't need to include rev comp.\n
  -unknown_low-quality_r1.fq\n
  -unknown_low-quality_r2.fq\n
  -idx_hop_r1.fq\n
  -idx_hop_r2.fq\n

-Open all of the files created above for writing.\n
-Within those, open the the 4 input files for reading.\n
-Create a for loop to simultaneously add the first record of r1, r2, r3, and r4 file to empty r1list, r2list, r3list, and r4list, respectively.\n

Examining index 1 currently in memory in the r2list:\n
-For r2list, check the index seq (pos 1 of r2list) and determine whether index is in dictionary.\n
  -If index is not in dictionary, it is unknown (Or could have low quality score, but don't need to check).\n
    -Call function str_add(r1list, r4list).\n
      -write cor_r1list to unknown_low-quality_r1.fq.\n
      -write cor_r4list to unknown_low-quality_r2.fq.\n
      -Wipe all lists and move to next record for all 4 input files\n
  -If index is in dictionary, determine quality score at each nucleotide in index 1 read and if all nucleotides are above the quality score coverage cut off, move to r3 file (index 2) following the code in the next block.

Examining index 2 currently in memory in the r3list:\n
-Should only have to run this block if index 1 was present in dictionary and had passing quality score for each nt in index 1.\n
-Call function convert_phred(letter) to determine if qscore of each nt in index 2 seq stored in the r3list is above coverage cutoff.\n
-If above coverage cutoff:\n
  -For r3list, check whether index seq (pos 1 of r3list) is the reverse compliment of the current index 1 seq in pos 1 of r2list. Call function rev_comp(sequence) to determine rev comp.\n
  -If index 2 is rev comp of index 1:\n
    -Call function str_add(r1list, r4list).\n
      -write cor_r1list to idx1seq_r1.fq.\n
      -write cor_r4list to idx1seq_r2.fq.\n
      -Wipe all lists and move to next record for all 4 input files\n
  -If index 2 not rev comp of index 1, check rev comp against the dict.\n  
    -If not in dict:\n
      -Call function str_add(r1list, r4list).\n
        -write cor_r1list to unknown_low-quality_r1.fq.\n
        -write cor_r4list to unknown_low-quality_r2.fq.\n
        -Wipe all lists and move to next record for all 4 input files\n
    -If in dict:\n
      -Call function str_add(r1list, r4list).\n
        -write cor_r1list to idx_hop_r1.fq.\n
        -write cor_r4list to idx_hop_r2.fq.\n
        -Wipe all lists and move to next record for all 4 input files\n
-If a single nt index 2 is below coverage cutoff:\n
  -Call function str_add(r1list, r4list).\n
    -write cor_r1list to unknown_low-quality_r1.fq.\n
    -write cor_r4list to unknown_low-quality_r2.fq.\n
    -Wipe all lists and move to next record for all 4 input files\n

-Open all of the newly created read 1 (or read 2) output files\n
-Use a for loop to iterate over each read in each file and count the number of reads. Or count the number of total lines and divide by 4.\n
-print the file name and number of reads calculated above for each read 1 output file.\n

-To determine each possible pair of indexes: The dual matched have already been determined, just need to determine all types of possible swapped indexes.\n
-Initialize a new empty dictionary\n
  -For just the read 1 idx_hop_r1.fq, add the index string appended to the header to the key of the dictionary along with the number of times this index string has been seen as the value.\n
  -Count through this file, continually adding occurrences of each string to dictionary.\n
  -Print out the dictionary.\n
  -Note, I will not do this for the unknown/low quality output files as the index strings of those reads have qscores that are too low.

FUNCTION for determining reverse compliment\n
-Function will take in a sequence and return its reverse compliment\n
def rev_comp(sequence):\n
return reverse_seq\n
Test parameter: "ATG"\n
Test output: "CAT"

FUNCTION for calculating qscore\n
-Function will calculate qscore from ASCI character\n
def convert_phred(letter):\n
return ord(letter) - 33\n
Test parameter: !\n
Test output: 0

FUNCTION for adding index string to header of biological reads\n
-Function will take the two biological read records currently being looked at and stored in memory and append the index 1 sequence as appears and the rev comp of index 2 seq to header of each biological read separated by a "-" (i.e. AAAAAAAA-AAAAAAAA)\n
def str_add(r1list, r4list):\n
return cor_r1list, cor_r4list\n
Test parameter:\n
"@K00337:83:HJKJNBBXX:8:1101:1347:1191 1:N:0:1,\n
GNGGTCTTCTACCTTTCTCTTCTTTTTTGGAGGAGTAGAATGTTGAGAGTCAGCAGTAGCCTCATCATCACTAGATGGCATTTCTTCTGAGCAAAACAGGT,\n
+,\n
A#A-<FJJJ<JJJJJJJJJJJJJJJJJFJJJJFFJJFJJJAJJJJ-AJJJJJJJFFJJJJJJFFA-7<AJJJFFAJJJJJF<F--JJJJJJF-A-F7JJJJ\n
"\n
Test output:\n
written in r1 file\n
"@K00337:83:HJKJNBBXX:8:1101:1347:1191 1:N:0:1_GTAGCGTA-GTAGCGTA,\n
GNGGTCTTCTACCTTTCTCTTCTTTTTTGGAGGAGTAGAATGTTGAGAGTCAGCAGTAGCCTCATCATCACTAGATGGCATTTCTTCTGAGCAAAACAGGT,\n
+,\n
A#A-<FJJJ<JJJJJJJJJJJJJJJJJFJJJJFFJJFJJJAJJJJ-AJJJJJJJFFJJJJJJFFA-7<AJJJFFAJJJJJF<F--JJJJJJF-A-F7JJJJ"\n
written in r2 file\n
"@K00337:83:HJKJNBBXX:8:1101:1347:1191 4:N:0:1_GTAGCGTA-GTAGCGTA,\n
NAAATGCCATCTAGTGATGATGAGGCTACTGCTGACTCTCAACATTCTACTCCTCCAAAAAAGAAGAGAAAGATTCCAACCCCCAGAACCGATGACCGGCA,\n
+,\n
#AAFFJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJFJJJJJJJJJJJFJJFJJJFJFJAJAF<7AJF<J--7AA7<FJ----7A-F-77-77------"
