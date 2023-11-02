#!/usr/bin/env python

# samtools view -b /projects/bgmp/shared/deduper/Dataset1.sam > Dataset1.bam 
# samtools sort Dataset1.bam -o Dataset1.sorted.bam
# samtools view -h Dataset1.sorted.bam > Dataset1.sorted.sam

#./deduper.py -u STL96.txt -f unittests/input.sam -o test.sam

import bioinfo
import re

#set up argparse with help statement
import argparse
def get_args():
    parser = argparse.ArgumentParser(description="deduper.py removes PCR duplicates. It takes two files: a file with UMIs and a SAM file you would like to deduplicate. It outputs two files: one deduplicated SAM file, and one SAM file containing all duplicates. The limitations of this function include: 1. Does not account for hard clipping 2. Does not use new CIGAR string encoding of '=', 'x'. Only uses 'M' 3. Ignores N's (splicing) for plus strand, so some duplicates are included in the final deduplicated output.")
    parser.add_argument("-u", "--umi", help="Specify the file name containing UMI list", required=True, type=str)
    parser.add_argument("-f", "--file", help="Specify the file SAM file you want to deduplicate", required=True, type=str)
    parser.add_argument("-o", "--outfile", help="Specify the name of the deduplicated output file", required=True, type=str)
    return parser.parse_args()
args = get_args()

#1a. Function for creating umi set
def umi_reference_set(umifile):
    umi_set = set()
    with open(umifile, "r") as fh:
        for line in fh:
            line = line.strip()
            umi_nucleotides = line.split()
            if bioinfo.validate_base_seq(line) == False:
                raise Exception("Unexpected character in your UMIs. UMI must contain only nucleotides ACTGN.")
            elif line in umi_set:
                raise Exception("Your UMIs are not all unique")
            else:
                umi_set.add(line)
    return umi_set

#1b. Make umi reference set from the file
umis = umi_reference_set(args.umi)

#setup splitting function for parsing lines in the sam file
def splitit(splitline):
    '''Takes a split line and returns a list of raw position, umi, strand, and cigar string.'''
    strand = ""
    #get umi
    qname = splitline[0]
    qname = qname.split(":")
    umi = qname[-1]
        
    #get strandedness 
    flag = int(splitline[1])
    if ((flag & 16)) == 16:
        strand = "minus"
    else:
        strand = "plus"

    #get the raw position 
    rawposition = splitline[3]

    #get the cigar string
    cigar = splitline[5]
    return rawposition, umi, strand, cigar

def adjust_plus(cigar, rawposition):
    number = ""
    for a in cigar:
        if a.isnumeric() == True:
            number = number + str(a)
        else:
            if a == "S": # we only have to look at first soft-clip event for plus strand reads
                adjustedposition = int(rawposition) - int(number)
                break #otherwise, if you have another softclip at the end, it will be wrong
            else:
                adjustedposition = int(rawposition)
                break
    return adjustedposition

def cigar_cutter(cigar: str):
    operator = []
    number = []
    split = re.split('(\d+)', cigar) #split into a list by occurence of a digit
     
    for a in split:
        if re.fullmatch(r"\d+", a):
            number.append(a)
        else:
            operator.append(a)
             
    operator = operator[1:] #splits into empty space at the beginning of the string, this removes the empty string.
    return number, operator

def adjust_minus(number: list, operator: list, rawposition: int):
    '''Retrieve the 5' start position of the minus strand, given its leftmost "raw" position. Accounts for M, D, N, S. Ignores I.'''
    adjustedposition = int(rawposition)
    if operator[0] == "S":
        adjustedoperators = operator[1:] #remove first soft clip
        adjustednumber = number[1:] #remove first soft clip
        for i,a in enumerate(adjustedoperators):
            if a == "M" or a == "D" or a == "N" or a == "S": #last S is in case there's soft clipping on end
                adjustedposition += int(adjustednumber[i])
            else: #a == I (which does not have an insertion to the reference), or any other cigar string which we do not account for
                adjustedposition = adjustedposition
    else:
         for i,a in enumerate(operator):
            if a == "M" or a == "D" or a == "N" or a == "S": #last S is in case there's soft clipping on end
                adjustedposition += int(number[i])
            else: #a == I (which does not have an insertion to the reference), or any other cigar string which we do not account for
                adjustedposition = adjustedposition
    return adjustedposition-1

#set up empty dictionary

#unique_reads is a set that contains strings of (rawposition:umi:strand)
unique_reads = set()
chromosome = ""

#Report
numberheaderlines = 0
numberuniquereads = 0
numberwrongumis = 0
removeddups = 0

#File names
samfile = open(args.file, "r")
outfile = open(args.outfile, "w")
duplicatefile = args.outfile + "_duplicates"
duplicateoutfile = open(duplicatefile, "w")
report = open("Deduper_Report.txt", "w")

for line in samfile:
    if line.startswith("@"):
        outfile.write(line)
        duplicateoutfile.write(line)
        numberheaderlines += 1
    else:
        splitline = line.split("\t")
        current_chromosome = splitline[2]
        #check if chromosome is the same
        if current_chromosome != chromosome:
            #empty the set at each new chromosome
            report.write(f"Chromosome {chromosome}, Unique reads: {len(unique_reads)}\n")
            unique_reads = set()
            chromosome = current_chromosome
            print(f"Now deduplicating reads on chromosome {chromosome}")
            rawposition, umi, strand, cigar = splitit(splitline)
            if umi not in umis:
                numberwrongumis += 1

            else:
                if strand == "plus":
                    adjustedposition = adjust_plus(cigar, rawposition)

                else: #minus strand
                    number, operator = cigar_cutter(cigar)
                    adjustedposition = adjust_minus(number, operator, rawposition)

                read_ID = str(adjustedposition) + ":" + umi + ":" + strand
                if read_ID in unique_reads:
                    duplicateoutfile.write(line)
                    removeddups +=1
                else:
                    outfile.write(line)
                    unique_reads.add(read_ID)
                    numberuniquereads += 1
        else:
            rawposition, umi, strand, cigar = splitit(splitline)
            if umi not in umis:
                numberwrongumis += 1
            else:
                if strand == "plus":
                    adjustedposition = adjust_plus(cigar, rawposition)

                else:
                    number, operator = cigar_cutter(cigar)
                    adjustedposition = adjust_minus(number, operator, rawposition)

                read_ID = str(adjustedposition) + ":" + umi + ":" + strand
                if read_ID in unique_reads:
                    duplicateoutfile.write(line)
                    removeddups +=1
                else:
                    outfile.write(line)
                    unique_reads.add(read_ID)
                    numberuniquereads += 1
                    # print(unique_reads)

report.write(f"Chromosome {chromosome}, Unique reads: {len(unique_reads)}\n")
report.write(f"Number of header lines = {numberheaderlines}\nNumber of unique reads = {numberuniquereads}\n")
report.write(f"Number of wrong umis = {numberwrongumis}\nNumber of duplicates = {removeddups}\n")

samfile.close()
outfile.close()
duplicateoutfile.close()
report.close()


