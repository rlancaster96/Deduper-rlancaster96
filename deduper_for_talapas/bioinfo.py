# #!/usr/bin/env python < Uncomment this to run assert statements

# Author: Ruben Lancaster

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions written during the Bioinformatics and Genomics Program coursework.'''

__version__ = "0.6.2"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = set("ATGCNatcgn")
RNA_bases = set("AUGCNaucgn")
revdict = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}

#1
def calc_median(sortedlist: list) -> int | float:
    '''Takes a sorted one-dimensional list and returns the median.'''
    n = len(sortedlist) 
    if n%2==0: #if the list is even:
        p1 = int(n/2)
        p2 = int((n/2)-1)
        median = (sortedlist[p1] + sortedlist[p2])/2
        return median
    else: # if the list is odd:
        p3 = int(n/2)
        median = sortedlist[p3]
        return median

#2
def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score. Uses a -33 phred convert.'''
    phred = int(ord(letter))-33  #Convert the letter to Unicode, convert to an integer, then subtract -33 to get the phred score 
    return phred

#3
def qual_score(phred_score: str) -> float:
    '''This function takes a phred score (as a string) and returns the average quality score.'''
    phred_tup = tuple(phred_score) # convert phred_score string to tuple.
    total = 0 #empty total
    for score in phred_tup:
        numericscore = convert_phred(score) # convert to numeric quality score
        total += numericscore # add score to running total
    return total/(len(phred_score))

#4
def validate_base_seq(seq,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNA_bases if RNAflag else DNA_bases)

#5
def gc_content(DNA):
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(DNA), "String contains invalid characters - are you sure you used a DNA sequence?"
    DNA = DNA.upper()
    return (DNA.count("G")+DNA.count("C"))/len(DNA)

#6
def oneline_fasta(readfile, writefile):
    '''Converts a fasta file with multiple sequence lines per record to one sequence line per record.'''
    with open(readfile, "r") as readfh, open(writefile, "w") as writefh:
        seq = ""
        while True:
            line = readfh.readline().strip()
            if not line:
                break
            if line.startswith(">"):
                if seq != "":
                    writefh.write(f"{seq}\n")
                seq = ""
                writefh.write(f"{line}\n")
            else:
                seq += line
        writefh.write(seq)
        return

#7 
def rev_comp(seq: str) -> str:
    '''Takes a dictionary called revdict (containing nucleotide complements) and a string. Returns reverse complement string. Revdict variable defined in bioinfo module.'''
    revlist = (list(reversed(seq)))
    revcomp = []
    for nuc in revlist:
        revcomp.append(revdict.get(nuc))
    revcomp = ''.join(revcomp)
    return revcomp

#8
def makedict(index_names: dict) -> dict:
    '''Takes a dictionary with index strings as values and returns a new dictionary with keys = index string, values = reverse compliment of index string.'''
    index_dict={}
    for index in index_names:
        indexrev = rev_comp(index_names.get(index))
        index_dict[index_names.get(index)] = indexrev
    return index_dict


if __name__ == "__main__":
    # Tests for calc_median
    assert calc_median([0,1,2,3,4]) == 2, "Wrong median for odd list [0,1,2,3,4]"
    assert calc_median([0,1,2,3,4,5]) == 2.5, "Wrong median for even list [0,1,2,3,4,5]"
    assert calc_median([-1,-3,-5]) == -3, "Wrong median for list with negative numbers [-1,-3,-5]"
    assert calc_median([0,0,0,0,0,0]) == 0, "Wrong median for list of the same value"
    print("Your calc_median function is working 1/5")

if __name__ == "__main__":
    #Tests for convert_phred (From Leslie)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working 2/5")


if __name__ == "__main__":
    # Tests for qual_score
    assert qual_score("@") == 31, "Your qual_score function does not convert phred scores correctly"
    assert qual_score("IIIIIIIIII") == 40, "Wrong average for string IIIIIIIIII"
    assert qual_score("$$$@@@") == 17, "Wrong average for string $$$@@@"
    print("Your qual_score function is working 3/5")

if __name__ == "__main__":
    # Tests for validate_base_seq
    assert validate_base_seq("AGCTGATCGA") == True, "validate_base_seq breaks for DNA"
    assert validate_base_seq("GACGUACGUAGCGCU", True) == True, "validate_base_seq breaks for RNA"
    assert validate_base_seq("189238403") == False, "validate_base_seq breaks for numbers"
    assert validate_base_seq("ACTCTACTATCAP") == False, "validate_base_seq breaks for other letters"
    assert validate_base_seq("ACTCTANNNTAGTGA") == True, "validate_base_seq breaks for letter N"
    assert validate_base_seq("agcuguauguau", True) == True, "validate_base_seq breaks for lowercase"
    print("Your validate_base_seq function is working 4/5")

if __name__ == "__main__":
    # Tests for gc_content
    assert gc_content("GGGGGGGGGGG") == 1, "gc_content breaks for G"
    assert gc_content("CCC") == 1, "gc_content breaks for C"
    assert gc_content("cccc") == 1, "gc_content breaks for lowercase"
    assert gc_content("AGACAGAC") == 0.5, "gc_content does not calculate correct average"
    print("Your gc_content function is working 5/5")
