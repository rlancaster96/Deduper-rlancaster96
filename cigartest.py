#!/usr/bin/env python

# cigar = "71M4S"

# def adjust_minus(rawposition, cigar):
#     cigardict = {}
#     allnumber = ""
#     for a in cigar:
#         if a.isnumeric() == True:
#             number = str(a)
#             allnumber = allnumber + number
#         else:
#             operator = a
#             cigardict[a] = allnumber
#             allnumber = ""
#     adjustedposition = int(rawposition) - int(cigardict['S'])
#     return adjustedposition

# position = adjust_minus(100, "20M2S")

# print(position)

import re

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
    adjustedposition = rawposition
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
    return adjustedposition

rawposition = 100
number, operator = cigar_cutter("50S10M10S")
adjustedposition = adjust_minus(number, operator, rawposition)

print(number)
print(operator)
print("raw: ", rawposition, "adjusted: ", adjustedposition)








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




# OLD MINUS 

# def adjust_minus(cigar, rawposition):
#     cigardict = {}
#     allnumber = ""
#     for a in cigar:
#         if a.isnumeric() == True:
#             number = str(a)
#             allnumber = allnumber + number
#         else:
#             operator = a
#             cigardict[a] = allnumber
#             allnumber = ""
#         #allows me to overwrite the first S, gets the last S which is the softclip closer to the 5' end.
#     adjustedposition = int(rawposition) - int(cigardict['S'])
#     return adjustedposition
