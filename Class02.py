# -*- coding: utf-8 -*-
"""
Created on Tue Jan 20 10:44:04 2015

@author: Marcola
"""


"""
With the skills you have now learned in the first part of this script, try this
exercise:
*** Sequence Manipulation Exercise ***
- Create a new Python script (text file)
- At the beginning of the script, define a DNA sequence (taken from 
https://github.com/jembrown/CompPhylo_Spr2015/blob/master/CodingSeq.txt)
- Print the length of the sequence to the screen along with text explaining 
the value
- Create and store the RNA equivalent of the sequence, then print to screen.
- Create and store the reverse complement of your sequence, then print to 
screen.
- Extract the bases corresponding to the 13rd and 14th codons from the 
sequence, then print them to the screen.
- Create a function to translate the nucleotide sequence to amino acids 
using the vertebrate mitochondrial genetic code (available from 
https://github.com/jembrown/CompPhylo_Spr2015/blob/master/VertMitTransTable.txt).
- Translate the sequence and print it to the screen.
- Be sure you've added comments to explain what this script is and what the 
different bits of code mean.
- Save this script as "seqManip.py" and commit it to your class GitHub repo.
"""

# defining the DNA sequence as a variable
dnaSeq = "aaaagctatcgggcccataccccaaacatgttggttaaaccccttcctttgctaattaatccttacgctatctccatcattatctccagcttagccctgggaactattactaccctatcaagctaccattgaatgttagcctgaatcggccttgaaattaacactctagcaattattcctctaataactaaaacacctcaccctcgagcaattgaagccgcaactaaatacttcttaacacaagcagcagcatctgccttaattctatttgcaagcacaatgaatgcttgactactaggagaatgagccattaatacccacattagttatattccatctatcctcctctccatcgccctagcgataaaactgggaattgccccctttcacttctgacttcctgaagtcctacaaggattaaccttacaaaccgggttaatcttatcaacatgacaaaaaatcgccccaatagttttacttattcaactatcccaatctgtagaccttaatctaatattattcctcggcttactttctacagttattggcggatgaggaggtattaaccaaacccaaattcgtaaagtcctagcattttcatcaatcgcccacctaggc"
print (dnaSeq)
#Defining how long is our DNA seq
print (len(dnaSeq))
dnaSeqLen = 618
print (dnaSeqLen)
#Replace "t" for "u" to qet the equivalent RNA sequence
rnaSeq = dnaSeq.replace("t","u")
print (rnaSeq)
#Creating the reverse complement of the DNA sequence
rev01=dnaSeq.replace("c","G") 
print (rev01)
rev02=rev01.replace("a","T")
print (rev02)
rev03=rev02.replace("t","A")
print (rev03)
rev04=rev03.replace("g","C")
print rev04
#Reversing the DNA sequence using the 'formula' [begin:end:step]
reverseDNA=rev04.lower()[::-1]
print (reverseDNA)
#- Extract the bases corresponding to the 13rd and 14th codons from the 
# sequence, then print them to the screen.
cod13 = dnaSeq[36:39]
print (cod13)
cod14 = dnaSeq[39:42]
print (cod14)
#Splitting the sequence into parts of 3 bases in order to get specific codons 
myCodons=[dnaSeq[x:x+3] for x in range(0,len(dnaSeq),3)]
for codons in myCodons:
    print codons
print myCodons
# The function 'codonMatch'  splits any dna sequence into codons, retrieving a 
# list with the right amino acid for each codon.
# The input dna sequence must be a string.
def codonMatch(dna):
    myCodons=[dna[x:x+3] for x in range(0,len(dna),3)] # splitting the codons
    for codon in myCodons: # Assigning each codon to the respective amino acid
        if codon in ['ttt','ttc']: 
            print "Phe",
        elif codon in ['tta', 'ttg', 'ctt', 'cta', 'ctc', 'ctg']:
            print "Leu",
        elif codon in ['att', 'atc']:
            print "Ile",
        elif codon in ['ata', 'atg']:
            print "Met",
        elif codon in ['gtt','gtc', 'gta', 'gtg']:
            print "Val",
        elif codon in ['tct', 'tcc', 'tca', 'tcg']:
            print "Ser",
        elif codon in ['cct', 'ccc', 'cca', 'ccg']:
            print "Pro",
        elif codon in ['act' , 'acc', 'aca' , 'acg']:
            print "Thr",
        elif codon in ['gct' , 'gcc', 'gca' , 'gcg']:
            print "Ala",
        elif codon in ['tat' , 'tac']:
            print "Tyr",
        elif codon in ['taa' , 'tag']:
            print "Ter",
        elif codon in ['cat' , 'cac']:
            print "His",
        elif codon in ['caa' , 'cag']:
            print "Gln",
        elif codon in ['aat' , 'aac']:
            print "Asn",
        elif codon in ['aaa' , 'aag']:
            print "Lys",
        elif codon in ['gat' , 'gac']:
            print "Asp",
        elif codon in ['gaa' , 'gag']:
            print "Glu",
        elif codon in ['tgt' , 'tgc']:
            print "Cys",
        elif codon in ['tga' , 'tgg']:
            print "Trp",
        elif codon in ['cgt' , 'cgc' , 'cga' , 'cgg']:
            print "Arg",
        elif codon in ['agt' , 'agc']:
            print "Ser",
        elif codon in ['aga' , 'agg']:
            print "Ter",
        elif codon in ['ggt' , 'ggc' , 'gga' , 'ggg']:
            print "Gly",
        else:
            print "sp"
            
codonMatch(dnaSeq)
    

        