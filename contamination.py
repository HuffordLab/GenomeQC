#!/usr/bin/env python
# coding: utf-8

"""
This script blasts the input genome assembly to the NCBI Univec database and generates a plot showing the number of scaffolds/contigs with blast hits

Usage: python contamination.py <genome assembly file (fasta format)> <output file name> <email address> 

"""

from __future__ import division
import time
import traceback
import sys
import os
import argparse
from argparse import RawTextHelpFormatter
import numpy as np
from Bio import SeqIO
import pandas as pd
import glob
import plotly.offline as offline
import plotly.graph_objs as go
from email.mime.text import MIMEText
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from smtplib import SMTP
import smtplib
import re
import subprocess
from Bio.Blast.Applications import NcbiblastnCommandline
from docopt import docopt
from collections import defaultdict
from os.path import basename, isfile, join, dirname, abspath
from sys import path
path.append(dirname(dirname(abspath(__file__))))
import lib.BtLog as BtLog
import lib.BtIO as BtIO




output_name = sys.argv[2] #name of the busco output directory
megablast_output = output_name.replace("_busco", "_megablast") #name of the megablast output file

#run megablast on the input genome assembly
cline = NcbiblastnCommandline(task="megablast", query=sys.argv[1], db="/home/nancy/assembly_app/blobtools/blobtools-master/Univec_modified.fasta", outfmt=6, max_target_seqs=1, max_hsps=1, num_threads=4, evalue=1e-25, out=megablast_output)

cline()

file_s = os.path.getsize(megablast_output)

#checks if the blast output has any hits or not

#if the blast output file is empty then the code generates barplot showing all the scaffolds/contigs with no blast hits 
if file_s == 0:
   records = list(SeqIO.parse(sys.argv[1], "fasta"))
   number_of_scaffolds = len(records)
   contigs_with_nohit = number_of_scaffolds

   d = {0: ['other_contigs']*contigs_with_nohit, 1: ['no_blast_hits']*contigs_with_nohit}
   
   #creates the dataframe
   df = pd.DataFrame(data=d)
   frames = [df]
   result = pd.concat(frames)
   result.columns = ['contig_name','species']
   df3 = result.groupby(['species']).size().reset_index(name='counts')
   df3.columns = ['contig_name', 'species']
   df4 = df3.sort_values(by=['species'], ascending=[True])
   xs = df4['contig_name']
   ys = df4['species']
   
   #plots the result of megablast
   data = [go.Bar(
            x=xs,
            y=ys,
             marker={

        'color': ys,

        'colorscale': 'Viridis'

      }
   )]


   layout = go.Layout(
   title='Contamination check against vector/adapator sequences',
   xaxis= dict(
        title=''
     ),
     yaxis=dict(
        title='sequence count'
       )
    
     )
   
   #name of the output plot file
   contamination_plot_name = output_name.replace("_busco", "contamination_plot.html")
   fig = go.Figure(data=data, layout=layout)
   offline.plot(fig, filename=contamination_plot_name, image='png')


   recipients = [str(sys.argv[3])] #email address provided by the user
   emaillist = [elem.strip().split(',') for elem in recipients]
   msg = MIMEMultipart()
   msg['Subject'] = "Here's your contamination plot"
   msg['From'] = 'assemblyshine@gmail.com'
   msg['Reply-to'] = 'assemblyshine@gmail.com'

   msg.preamble = 'Multipart massage.\n'
   
   #body of the email
   part = MIMEText("Please find attached the contamination plot. No blast hits were found in your fasta file. For contamination analysis, the assembled genome sequences are matched against the NCBI UniVec database. UniVec (https://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/) is a database that can be used to quickly identify segments within nucleic acid sequences which may be of vector origin (vector contamination). In addition to vector sequences, UniVec also contains sequences for those adapters, linkers, and primers commonly used in the process of cloning cDNA or genomic DNA. This enables contamination with these oligonucleotide sequences to be found during the vector screen.")
   msg.attach(part)

   part = MIMEApplication(open(contamination_plot_name,"rb").read())
   
   #attaches the bar chart html file to the email
   part.add_header('Content-Disposition', 'attachment', filename=contamination_plot_name)
   msg.attach(part)
   
   server = smtplib.SMTP("smtp.gmail.com:587")
   server.ehlo()
   server.starttls()
   server.login("assemblyshine@gmail.com", "rshiny@2020")


   server.sendmail(msg['From'], emaillist , msg.as_string())


#if blast hits are found in the input genome assembly, the code generates a barplot showing the number of scaffolds/contigs in the assembly with blast hits
else:
  #This code adds taxonID information to the sequence IDs (of the subject sequences) of the megablast output file. More information on this can be found here: https://blobtools.readme.io/docs/taxify 

  def mapping():
      out_f, hit_f, map_f, taxid_d = None, None, None, {}
      hit_f = megablast_output  #hit file: BLAST similarity search result (TSV format)
      map_f = "/home/nancy/assembly_app/blobtools/blobtools-master/taxon_n"  #mapping file (TSV format), in which one column lists a sequence ID (of a subject) and another the NCBI TaxID
      map_col_sseqid = "0" #column of mapping file containing sequence IDs (of the subject)
      map_col_taxid = "2"  #column of mapping file containing the TaxID of the subject
      hit_col_qseqid = "0" #column of the hit file containing query ID
      hit_col_sseqid = "1" #column of the hit file containing subject ID
      hit_col_score = "11" #column of the hit file containing (bit)score

    
      try:
         hit_col_qseqid = int(hit_col_qseqid)
         hit_col_sseqid = int(hit_col_sseqid)
         hit_col_score = int(hit_col_score)
      except ValueError:
          BtLog.error('41' % ("--hit_column_qseqid, --hit_column_sseqid and --hit_column_score"))

   
      if map_f:
          if map_col_sseqid and map_col_taxid:
              try:
                 map_col_sseqid = int(map_col_sseqid)
                 map_col_taxid = int(map_col_taxid)
              except ValueError:
                  BtLog.error('44')
              print BtLog.status_d['1'] % ("Mapping file", map_f)
              taxid_d = BtIO.parseDict(map_f, map_col_sseqid, map_col_taxid)
              out_f = BtIO.getOutFile("taxified",hit_f,"out")
          else:
              BtLog.error('44')
      else:
          BtLog.error('41')

      output = []
      print BtLog.status_d['1'] % ("similarity search result", hit_f)
      with open(hit_f) as fh:
          for idx, line in enumerate(fh):
              col = line.rstrip("\n").split()
              qseqid = col[hit_col_qseqid]
              sseqid = col[hit_col_sseqid]
              score = col[hit_col_score]
              tax_id = None
              if sseqid not in taxid_d:
                 BtLog.warn_d['12'] % (sseqid, map_f)
              tax_id = taxid_d.get(sseqid, "N/A")
              output.append("%s\t%s\t%s\t%s" % (qseqid, tax_id, score, sseqid))
      if output:
         with open(out_f, "w") as fh:
             print BtLog.status_d['24'] % out_f
             fh.write("\n".join(output) + "\n")

  if __name__ == '__main__':
         mapping()


  taxify_output = megablast_output + ".taxified.out"
  
  #creates the dataframe
  df = pd.read_table(taxify_output, header=None)
  df1 = df.drop(columns=[2,3])

  records = list(SeqIO.parse(sys.argv[1], "fasta"))
  number_of_scaffolds = len(records)
  nrow = len(df1.index)
  contigs_with_nohit = number_of_scaffolds - nrow
  
  
  d = {0: ['other_contigs']*contigs_with_nohit, 1: ['no_blast_hits']*contigs_with_nohit}
  df2 = pd.DataFrame(data=d)
  frames = [df1, df2]
  result = pd.concat(frames)
  result.columns = ['contig_name','species']
  df3 = result.groupby(['species']).size().reset_index(name='counts')
  df3.columns = ['contig_name', 'species']
  df4 = df3.sort_values(by=['species'], ascending=[True])
  xs = df4['contig_name']
  ys = df4['species']
  
  #plots the result of megablast
  data = [go.Bar(
            x=xs,
            y=ys,
             marker={

        'color': ys,

        'colorscale': 'Viridis'

        }
    )]


  layout = go.Layout(
  title='Contamination check against vector/adapator sequences',
  xaxis= dict(
        title=''
      ),
  yaxis=dict(
        title='sequence count'
        )
    )
  
  
  ##name of the output plot file
  contamination_plot_name = output_name.replace("_busco", "contamination_plot.html")
  fig = go.Figure(data=data, layout=layout)
  offline.plot(fig, filename=contamination_plot_name, image='png')


  recipients = [str(sys.argv[3])] #email address provided by the user
  emaillist = [elem.strip().split(',') for elem in recipients]
  msg = MIMEMultipart()
  msg['Subject'] = "Here's your contamination plot"
  msg['From'] = 'assemblyshine@gmail.com'
  msg['Reply-to'] = 'assemblyshine@gmail.com'

  msg.preamble = 'Multipart massage.\n'
  
  #body of the email
  part = MIMEText("Please find attached the contamination plot. For contamination analysis, the assembled genome sequences are matched against the NCBI UniVec database. UniVec (https://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/) is a database that can be used to quickly identify segments within nucleic acid sequences which may be of vector origin (vector contamination). In addition to vector sequences, UniVec also contains sequences for those adapters, linkers, and primers commonly used in the process of cloning cDNA or genomic DNA. This enables contamination with these oligonucleotide sequences to be found during the vector screen.")
  msg.attach(part)

  part = MIMEApplication(open(contamination_plot_name,"rb").read())
  
  #attaches the bar chart html file to the email
  part.add_header('Content-Disposition', 'attachment', filename=contamination_plot_name)
  msg.attach(part)


  server = smtplib.SMTP("smtp.gmail.com:587")
  server.ehlo()
  server.starttls()
  server.login("assemblyshine@gmail.com", "rshiny@2020")

  server.sendmail(msg['From'], emaillist , msg.as_string())               
