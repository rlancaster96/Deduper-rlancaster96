#!/usr/bin/env python

import bioinfo

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
    umi = qname[7]
        
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
            if a == "S": # we only have to look at soft-clipping for plus strand reads
                adjustedposition = int(rawposition) - int(number)
                break #otherwise, if you have another softclip at the end, it will be wrong
            else:
                adjustedposition = int(rawposition)
                break
    return adjustedposition

def adjust_minus(cigar, rawposition):
    cigardict = {}
    allnumber = ""
    for a in cigar:
        if a.isnumeric() == True:
            number = str(a)
            allnumber = allnumber + number
        else:
            operator = a
            cigardict[a] = allnumber
            allnumber = ""
    adjustedposition = int(rawposition) - int(cigardict['S'])
    return adjustedposition

#set up empty dictionary

#unique_reads is a set that contains strings of (rawposition:umi:strand)
unique_reads = set()
chromosome = ""


samfile = open(args.file, "r")
outfile = open(args.outfile, "w")
duplicatefile = args.outfile + "_duplicates"
duplicateoutfile = open(duplicatefile, "w")

for line in samfile:
    if line.startswith("@"):
        outfile.write(line)
        duplicateoutfile.write(line)
    else:
        splitline = line.split("\t")
        current_chromosome = splitline[2]
        #check if chromosome is the same
        if current_chromosome != chromosome:
            #empty the set at each new chromosome
            unique_reads = set()
            print("the set is reset")
            chromosome = current_chromosome
            print(f"The chromosome has been updated to {chromosome}")
            rawposition, umi, strand, cigar = splitit(splitline)
            if strand == "plus":
                adjustedposition = adjust_plus(cigar, rawposition)
            else: #minus strand
                adjustedposition = adjust_minus(cigar,rawposition)
            read_ID = str(adjustedposition) + ":" + umi + ":" + strand
            if read_ID in unique_reads:
                duplicateoutfile.write(line)
            else:
                outfile.write(line)
                unique_reads.add(read_ID)
            print(unique_reads)
        else:
            rawposition, umi, strand, cigar = splitit(splitline)
            if strand == "plus":
                adjustedposition = adjust_plus(cigar, rawposition)
            else:
                adjustedposition = adjust_minus(cigar,rawposition)
            read_ID = str(adjustedposition) + ":" + umi + ":" + strand
            if read_ID in unique_reads:
                duplicateoutfile.write(line)
            else:
                outfile.write(line)
                unique_reads.add(read_ID)
            
            unique_reads.add(read_ID)
            print(unique_reads)

samfile.close()
outfile.close()

