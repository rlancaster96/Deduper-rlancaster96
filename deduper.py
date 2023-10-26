#!/usr/bin/env python

#import bioinfo.py
# import bioinfo

#By chromosome: UMI | adjusted 5' position | strand |

#set up argparse 
import argparse
def get_args():
    parser = argparse.ArgumentParser(description="deduper description")
    parser.add_argument("-u", "--umi", help="Specify the file name containing UMI list", required=True, type=str)
    parser.add_argument("-f", "--file", help="Specify the file SAM file you want to deduplicate", required=True, type=str)
    parser.add_argument("-o", "--outfile", help="Specify the name of the deduplicated output file", required=True, type=str)
    return parser.parse_args()
args = get_args()

#Make function for creating umi set
def umi_reference_set(umifile):
    umi_set = set()
    with open(umifile, "r") as fh:
        for line in fh:
            line = line.strip()
            if line in umi_set:
                raise Exception("Your UMIs are not all unique")
            else:
                umi_set.add(line)
    return umi_set

#make umi reference set from the file
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


#set up empty dictionary

#unique_reads is a set that contains strings of (rawposition:umi:strand:cigar)
unique_reads = set()
chromosome = ""


samfile = open(args.file, "r")
outfile = open(args.outfile, "w")
for line in samfile:
    if line.startswith("@"):
        outfile.write(line)
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
            read_ID = rawposition + ":" + umi + ":" + strand + ":" + cigar
            unique_reads.add(read_ID)
            print(unique_reads)
        else:
            rawposition, umi, strand, cigar = splitit(splitline)
            read_ID = rawposition + ":" + umi + ":" + strand + ":" + cigar
            unique_reads.add(read_ID)
            print(unique_reads)

samfile.close()
outfile.close()

