#!/usr/bin/python

"""
This script calculates the NG values at various thresholds (1-100%) for the given input genome assembly.

Usage: python3/Python-3.6.4/python NG.py <genome assembly file (fasta format)> <output file name> <estimated genome size (in Mb)> 

"""

from Bio import SeqIO

import sys

import statistics

import numpy as np



inputfile = sys.argv[1] 
outputfile = sys.argv[2]
estimated_genome_size = float(sys.argv[3])


sys.stdout=open(outputfile,"w")

#reads all the sequences of the input file 
records = list(SeqIO.parse(inputfile, "fasta"))

#calculates the length of each sequence
len_seq = [len(rec) for rec in records]

#orders the sequences in the decreasing order of length
sorted_len = sorted(len_seq, reverse=True)


NX = []

for i in range(1,100):  #calculates for all thresholds 1-100%

    testSumNG = 0       

    NG = 0              

    NGcon = 0           

    for conNG in sorted_len:

        testSumNG += conNG   

        NGcon += 1

        if estimated_genome_size*i*1000000/100 < testSumNG: #checks when the sum of lengths of sequences becomes greater than the ith% of estimated genome size (since the input estimated genome size is in megabases therefore multiplying by 10^6 to convert to bp)

           NG = conNG   # sequence length which when added to the sum of lengths crosses the ith threshold of estimated genome size gives the NG value at that threshold

           NX.append(NG)

           break





for elem in NX:
        print (elem) 

           

