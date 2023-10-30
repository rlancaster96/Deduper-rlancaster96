#!/usr/bin/env python

cigar = "71M4S"

def adjust_minus(rawposition, cigar):
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

position = adjust_minus(100, "20M2S")

print(position)








# def adjust_plus(cigar, rawposition):
#     number = ""
#     for a in cigar:
#         if a.isnumeric() == True:
#             number = number + str(a)
#         else:
#             if a == "S": # we only have to look at soft-clipping for plus strand reads
#                 adjustedposition = rawposition - int(number)
#                 break #otherwise, if you have another softclip at the end, it will be wrong
#             else:
#                 adjustedposition = rawposition
#                 break
#     return adjustedposition

# adjustedposition = adjust_plus("71M50S", 100)
# print(adjustedposition)
