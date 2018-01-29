# FastaMutator
Given one fasta sequence, creates multifasta files containing mutant fasta sequences

## Description: 
Given one fasta sequence, creates multifasta files containing mutant fasta sequences
The script allows only point mutations. See EMBOSS::msbar (not R) for more complex mutation schemes.

## Explanations about filenames generated:
   e.g. __Fasta1_MR0.01_Nmut100.fasta__ means
__Fasta1__ original name 
__MR0.01__ mutation rate of 0.01 = 1% = (seqlength of 100 bases x 0.01 = 1 point mutation was introduced)
__Nmut100__ 100 mutated sequences. (the number at the end of each fasta header indicates the mutant id e.g. …._1, …_2

## Warning message: 
"Caution: Cannot create non-redundant mutants for Fasta1.txt and MR=0.01sequence= 622!!!"
This indicates that the algorithm could not find a unique mutant (the same mutants occur at least twice in this set)
Change "Ntries" to make more (lengthier) search of unique sequences.

## PARAMETERS (hard coded for now)
FILES=c("Fasta1.txt")
MutationRate=c(0.01,0.02,0.03,0.05,0.10)
Nmutants=100 # number of mutated fasta sequences to create
Ntries=10 # how many times to try searching for a new sequence of a variant?
