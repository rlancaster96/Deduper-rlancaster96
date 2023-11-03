#!/usr/bin/env python

#import modules
import bioinfo
import re
import gzip

#set up argparse with help statement
import argparse
def get_args():
    parser = argparse.ArgumentParser(description="deduper.py removes PCR duplicates. It takes two files: a file with UMIs and a sorted SAM file you would like to deduplicate. It outputs one SAM file containing only unique reads (non-PCR duplicates) NOT in compressed format. The limitations of this function include: 1. Does not account for hard clipping 2. Does not use new CIGAR string encoding of '=', 'x'. (Uses 'M') 3. Does not account for splicing (so there may be duplicates falsely included in the final deduplicated output).")
    parser.add_argument("-u", "--umi", help="Specify the file name containing UMI list. One UMI per line.", required=True, type=str)
    parser.add_argument("-f", "--file", help="Specify the SAM file you want to deduplicate. Takes .sam or .sam.gz", required=True, type=str)
    parser.add_argument("-o", "--outfile", help="Specify the name of the deduplicated output file", required=True, type=str)
    return parser.parse_args()
args = get_args()

#****************************** Functions ***********************************

def umi_reference_set(umifile: str):
    '''Takes a file with one UMI per line and returns a set of UMIs.'''
    umi_set = set()
    with open(umifile, "r") as fh:
        for line in fh:
            line = line.strip()
            if bioinfo.validate_base_seq(line) == False:
                raise Exception("Unexpected character in your UMIs. UMI must contain only nucleotides ACTGN.")
            elif line in umi_set:
                raise Exception("Your UMIs are not all unique")
            else:
                umi_set.add(line)
    return umi_set

def splitit(splitline: str):
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

def adjust_plus(cigar: str, rawposition):
    '''Takes a cigar string and raw position and adjusts the position. Only for plus strand.'''
    # This function looks at if the first operator is an S.
    # If it is, then it subtracts the digits before the S. If not, then no change is made.
        
    number = ""
    for a in cigar:
        if a.isnumeric() == True:
            number = number + str(a)
        else:
            if a == "S": # we only have to look at first soft-clip event for plus strand reads
                adjustedposition = int(rawposition) - int(number)
                break #break otherwise, if you have another softclip at the end, it will be wrong
            else:
                adjustedposition = int(rawposition)
                break
    return adjustedposition

def cigar_cutter(cigar: str):
    '''Split cigar string into a list of digits and a list of operators (including M, I, D, N, S)'''
    operator = []
    number = []
    split = re.split('(\d+)', cigar) #split into a list by occurence of a digit
     
    for a in split:
        if re.fullmatch(r"\d+", a):
            number.append(a) #creates a list of length n elements of each digit in the cigar string
        else:
            operator.append(a) #creates a list of length n of elements of each operator (letter) in the cigar string
             
    operator = operator[1:] #removes empty string at the beginning of the list
    return number, operator

def adjust_minus(number: list, operator: list, rawposition):
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
    return adjustedposition-1 #subtract one from position to account for counting 

def read_IDer(rawposition, umi, strand, cigar):
    '''takes elements of the line and returns a read ID with corrected position using adjustment and cigar cutting functions.'''
    if strand == "plus":
        adjustedposition = adjust_plus(cigar, rawposition)
    else: #minus strand
        number, operator = cigar_cutter(cigar)
        adjustedposition = adjust_minus(number, operator, rawposition)
    read_ID = str(adjustedposition) + ":" + umi + ":" + strand
    return read_ID


#******************************** Script *************************************

if __name__ == "__main__":
    #Make umi reference set from the file
    umis = umi_reference_set(args.umi)

    #initialize global variables 
    unique_reads = set() #unique_reads is a set that contains strings of (rawposition:umi:strand)
    chromosome = ""
    seen_before = set() #set of chromosomes already seen to keep track if we have a properly sorted sam file.

    #initialize report variables
    numberheaderlines = 0
    numberuniquereads = 0
    numberwrongumis = 0
    removeddups = 0

    #Open files
    samname = args.file
    samname = samname.split(".")
    if "sam" not in samname:
        raise Exception("Please provide a file in .sam or .sam.gz format.")
    if samname[-1] == "gz": #identify if sam file is in compressed format
        samfile = gzip.open(args.file, "rt")
    else:
        samfile = open(args.file, "r")


    outfile = open(args.outfile, "w")
    report = open("Deduper_Report.txt", "w")

    for line in samfile:
        if line.startswith("@"):
            outfile.write(line)

            numberheaderlines += 1
        else:
            splitline = line.split("\t")
            current_chromosome = splitline[2]

            #check if chromosome is the same
            if current_chromosome != chromosome:
                if current_chromosome in seen_before:
                    raise Exception("Your sam file has not been sorted correctly.")
                if chromosome != "": #don't print the first blank
                    report.write(f"Chromosome {chromosome}, Unique reads: {len(unique_reads)}\n")
                seen_before.add(current_chromosome)
                unique_reads = set() #empty the set at each new chromosome
                chromosome = current_chromosome #update the chromosome
                print(f"Now deduplicating reads on chromosome {chromosome}")
                rawposition, umi, strand, cigar = splitit(splitline)
                if umi not in umis:
                    numberwrongumis += 1
                else:
                    read_ID = read_IDer(rawposition, umi, strand, cigar)
                    if read_ID in unique_reads:
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
                    read_ID = read_IDer(rawposition, umi, strand, cigar)
                    if read_ID in unique_reads:
                        removeddups +=1
                    else:
                        outfile.write(line)
                        unique_reads.add(read_ID)
                        numberuniquereads += 1

    report.write(f"Chromosome {chromosome}, Unique reads: {len(unique_reads)}\n")
    report.write(f"Number of header lines = {numberheaderlines}\nNumber of unique reads = {numberuniquereads}\n")
    report.write(f"Number of wrong umis = {numberwrongumis}\nNumber of duplicates = {removeddups}\n")

    samfile.close()
    outfile.close()
    report.close()




