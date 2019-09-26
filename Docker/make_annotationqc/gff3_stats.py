#!/usr/bin/python

"""
This script calculates all the basic metrics for the given input genome annotation set.

Usage: python3/Python-3.6.4/python gff3_stats.py <genome annotation file (gff or gff3 format)> <output file name>  

"""


import sys
import statistics

from collections import Counter


inputfile = sys.argv[1]
outputfile = sys.argv[2]

sys.stdout=open(outputfile,"w")


num_gene_models = 0
gene_length_list = []
num_exons = 0
num_transcripts = 0
exon_length_list = []
gene_models_less200 = []
exons = []
token = []
for line in open(inputfile,'r'):
      line = line.strip()
      if len(line)>0:
          token = line.split('\t')
          if len(token) > 4:
              if token[2] == "gene":
                  num_gene_models+=1
                  gene_length = (int(token[4]) - int(token[3]))+1
                  gene_length_list.append(gene_length)
                  if gene_length < 200:
                     gene_models_less200.append(gene_length)
              if token[2] == "exon":
                  num_exons+=1
                  exon_length = (int(token[4]) - int(token[3]))+1
                  exon_length_list.append(exon_length)
                  column9_exon = token[8]
              if token[2] == "mRNA" or token[2] == "transcript":
                 num_transcripts+=1 
                 


print("Number of gene models:", num_gene_models)
print("Minumum gene length:", min(gene_length_list))
print("Maximum gene length:", max(gene_length_list))
average_gene_length = round(sum(gene_length_list)/len(gene_length_list),1)
print("Average gene length:", average_gene_length)
print("Total number of exons:", num_exons)
average_number_exons = round((num_exons/num_gene_models),1)
print("Average number of exons per gene model:", average_number_exons)
average_exon_length = round((sum(exon_length_list)/len(exon_length_list)),1)
print("Average exon length:", average_exon_length)
print("Number of transcripts:", num_transcripts)
average_number_transcripts = round((num_transcripts/num_gene_models),1)
print("Average number of transcripts per gene model:",average_number_transcripts)
print("Number of gene models less than 200 bp:", len(gene_models_less200))


