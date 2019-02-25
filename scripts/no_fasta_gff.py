#!/usr/bin/env python
# coding: utf-8

"""
This script plots and email the scores for the user selected reference genome assemblies and annotations

Usage: python no_fasta_gff.py <email address> <busco assembly output directory name> <busco annotation output directory name> <reference genome names>

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






ouput_name_assembly = sys.argv[2] #name of busco assembly directory
ouput_name_gff = sys.argv[3]      #name of busco annotation directory
ref_names = sys.argv[4].split(',') #name of the reference genomes


#plotting busco scores for reference genome assemblies

L = []
for i in ref_names:
    text_file = open(("/home/nancy/assembly_app/reference_genomes/" + i + "_buscometrics"), 'r')
    L.append(text_file.readlines())


myarray = np.asarray(L)
print (myarray)
df = pd.DataFrame(data=myarray)
df = df.replace('\n','', regex=True)
df1 = df.T
df1.columns = ref_names
busco_categories = ["C&S", "D", "F", "M"]
df1.index = busco_categories
t1 = np.asarray(df1.iloc[0])
t2 = np.asarray(df1.iloc[1])
t3 = np.asarray(df1.iloc[2])
t4 = np.asarray(df1.iloc[3])
colnames = list(df1)

trace1 = go.Bar(
    x=colnames,
    y=t1,
    name='C&S: complete and single copy'
)
trace2 = go.Bar(
    x=colnames,
    y=t2,
    name='D: complete and duplicate copy'
)

trace3 = go.Bar(
    x=colnames,
    y=t3,
    name='F: fragmented'
)
trace4 = go.Bar(
    x=colnames,
    y=t4,
    name='M: missing'
)

data = [trace1, trace2, trace3, trace4]
layout = go.Layout(
    barmode='stack',
    title='Busco Assembly Plot',
        xaxis=dict(
        title='Genomes',
    ),
    yaxis=dict(
        title='Percentage of busco genes',
    )
)

busco_plot_name =  ouput_name_assembly + "assembly_plot.html"
fig = go.Figure(data=data, layout=layout)
offline.plot(fig, filename=busco_plot_name, image='png')

recipients = [str(sys.argv[1])]
emaillist = [elem.strip().split(',') for elem in recipients]
msg = MIMEMultipart()
msg['Subject'] = "Here's your busco assembly plot"
msg['From'] = 'assemblyshine@gmail.com'
msg['Reply-to'] = 'assemblyshine@gmail.com'

msg.preamble = 'Multipart massage.\n'

part = MIMEText("Hi, please find the attached busco plot. BUSCO scores are used to assess genome assembly and annotation completeness with Benchmarking Universal Single-Copy Orthologs. A high quality genome assembly should have higher number of Complete and Single copy genes and lower number of Fragmented and Missing genes. To get more information on the these busco analysis, please visit this link: https://busco.ezlab.org/")
msg.attach(part)

part = MIMEApplication(open(busco_plot_name,"rb").read())
part.add_header('Content-Disposition', 'attachment', filename=busco_plot_name)


msg.attach(part)


server = smtplib.SMTP("smtp.gmail.com:587")
server.ehlo()
server.starttls()
server.login("assemblyshine@gmail.com", "rshiny@2020")
server.sendmail(msg['From'], emaillist , msg.as_string())



#plotting busco scores for reference genome annotations

L1 = []
for i in ref_names:
    text_file_gff = open(("/home/nancy/assembly_app/reference_genomes/" + i + "_gffbuscometrics"), 'r')
    L1.append(text_file_gff.readlines())

myarray_gff = np.asarray(L1)
print (myarray_gff)
df_gff = pd.DataFrame(data=myarray_gff)
df_gff = df_gff.replace('\n','', regex=True)
df1_gff = df_gff.T
df1_gff.columns = ref_names
busco_categories = ["C&S", "D", "F", "M"]
df1_gff.index = busco_categories
t1 = np.asarray(df1_gff.iloc[0])
t2 = np.asarray(df1_gff.iloc[1])
t3 = np.asarray(df1_gff.iloc[2])
t4 = np.asarray(df1_gff.iloc[3])
colnames = list(df1_gff)

trace1 = go.Bar(
    x=colnames,
    y=t1,
    name='C&S: complete and single copy'
)
trace2 = go.Bar(
    x=colnames,
    y=t2,
    name='D: complete and duplicate copy'
)

trace3 = go.Bar(
    x=colnames,
    y=t3,
    name='F: fragmented'
)
trace4 = go.Bar(
    x=colnames,
    y=t4,
    name='M: missing'
)

data_gff = [trace1, trace2, trace3, trace4]
layout = go.Layout(
    barmode='stack',
    title='Busco Annotation Plot',
        xaxis=dict(
        title='Genomes',
    ),
    yaxis=dict(
        title='Percentage of busco genes',
    )
)

buscogff_plot_name =  ouput_name_gff + "_plot.html"
fig = go.Figure(data=data_gff, layout=layout)
offline.plot(fig, filename=buscogff_plot_name, image='png')

recipients = [str(sys.argv[1])]
emaillist = [elem.strip().split(',') for elem in recipients]
msg = MIMEMultipart()
msg['Subject'] = "Here's your busco annotation plot"
msg['From'] = 'assemblyshine@gmail.com'
msg['Reply-to'] = 'assemblyshine@gmail.com'

msg.preamble = 'Multipart massage.\n'

part = MIMEText("Hi, please find the attached busco plot. BUSCO scores are used to assess genome assembly and annotation completeness with Benchmarking Universal Single-Copy Orthologs. A high quality genome annotations should have higher number of Complete and Single copy genes and lower number of Fragmented and Missing genes. To get more information on the these busco analysis, please visit this link: https://busco.ezlab.org/")
msg.attach(part)

part = MIMEApplication(open(buscogff_plot_name,"rb").read())
part.add_header('Content-Disposition', 'attachment', filename=buscogff_plot_name)


msg.attach(part)


server = smtplib.SMTP("smtp.gmail.com:587")
server.ehlo()
server.starttls()
server.login("assemblyshine@gmail.com", "rshiny@2020")
server.sendmail(msg['From'], emaillist , msg.as_string())

