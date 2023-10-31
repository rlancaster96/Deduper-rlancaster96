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

def cigar_cutter(cigar):
    operator = []
    number = []
    split = re.split('(\d+)', cigar) #split into a list of digits
     
    for a in split:
        if re.fullmatch(r"\d+", a):
            number.append(a)
        else:
            operator.append(a)
             
    operator = operator[1:]
    return number, operator

number, operator = cigar_cutter("2S40M2I40M2S")

print(number)
print(operator)
 








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
