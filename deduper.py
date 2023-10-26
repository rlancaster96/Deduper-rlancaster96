#!/usr/bin/env python

#import bioinfo.py
# import bioinfo

#By chromosome: UMI | adjusted 5' position | strand |

#set up argparse 
import argparse
def get_args():
    parser = argparse.ArgumentParser(description="deduper description")
    parser.add_argument("-U", "--umis", help="Specify the file name containing UMI list", required=True, type=str)
    parser.add_argument("-S", "--samfile", help="Specify the file SAM file you want to deduplicate", required=True, type=str)
    parser.add_argument("-O", "--outfile", help="Specify the file SAM file you want to deduplicate", required=True, type=str)
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
umis = umi_reference_set(args.umis)

samfile = open(args.samfile, "r")
outfile = open(args.outfile, "w")
for line in samfile:
    if line.startswith("@"):
        outfile.write(line)
    else:
        outfile.write("THIS LINE WILL BE ASSESSED FOR DEDUPER\n")

samfile.close()
outfile.close()

