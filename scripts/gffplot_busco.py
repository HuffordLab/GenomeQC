#!/usr/bin/env python
# coding: utf-8


"""
This script plots and email the scores for the user uploaded genome annotations and compare the scores to the user selected reference genomes

Usage: python gffplot_busco.py <email address> <name for output plot> <busco output directory name> <reference genome names>

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





L = []
ref_names = sys.argv[4].split(',') #name of the reference genomes
output_name = sys.argv[3]          #name of busco directory
user_n  = sys.argv[2]              #user provided name to label the plot 
name_busco_directory = "/opt/shiny-server/samples/sample-apps/assembly_statistics/run_" + output_name

summary_file = name_busco_directory + "/" + "short_summary_" + output_name + ".txt"

#reads the calculated busco scores from the summary file for the user uploaded genome annotations   
f=open(summary_file)
lines=f.readlines()
b = lines[7]
c = re.findall(r'[\d\.\d]+', b)
user_input = c[1:-1]


#reads the busco scores of the reference genomes
for i in ref_names:
    text_file = open(("/home/nancy/assembly_app/reference_genomes/" + i + "_gffbuscometrics"), 'r')
    L.append(text_file.readlines())



#creates the panda dataframe
myarray = np.asarray(L)
print (myarray)
df = pd.DataFrame(data=myarray)
df = df.replace('\n','', regex=True)
df1 = df.T
df1.columns = ref_names
df1 = df1.assign(user_name=user_input)
df2 = df1.rename(columns = {'user_name':user_n})
busco_categories = ["C&S", "D", "F", "M"]
df2.index = busco_categories

t1 = np.asarray(df2.iloc[0])
t2 = np.asarray(df2.iloc[1])
t3 = np.asarray(df2.iloc[2])
t4 = np.asarray(df2.iloc[3])
colnames = list(df2)

#plots the busco scores as stacked bar chart
trace1 = go.Bar(
    x=colnames,
    y=t1,
    name='C&S: complete & single copy'
)
trace2 = go.Bar(
    x=colnames,
    y=t2,
    name='D: complete & duplicate copy'
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
    title='Busco Plot',
        xaxis=dict(
        title='Genomes',

    ),
    yaxis=dict(
        title='Percentage of busco genes',

    )
)

#name of the output plot file
busco_plot_name = output_name + "_plot.html"
fig = go.Figure(data=data, layout=layout)
offline.plot(fig, filename=busco_plot_name, image='png')

recipients = [str(sys.argv[1])] #email address provided by the user
emaillist = [elem.strip().split(',') for elem in recipients]
msg = MIMEMultipart()
msg['Subject'] = "Here's your annotation busco plot"
msg['From'] = 'assemblyshine@gmail.com'
msg['Reply-to'] = 'assemblyshine@gmail.com'

msg.preamble = 'Multipart massage.\n'

#body of the email
part = MIMEText("Hi, please find the attached busco plot. BUSCO scores are used to assess genome assembly and annotation completeness with Benchmarking Universal Single-Copy Orthologs. A high quality genome annotations should have higher number of Complete and Single copy genes and lower number of Fragmented and Missing genes. To get more information on the these busco analysis, please visit this link: https://busco.ezlab.org/")
msg.attach(part)

part = MIMEApplication(open(busco_plot_name,"rb").read())

#attaches the stacked bar chart html file to the email
part.add_header('Content-Disposition', 'attachment', filename=busco_plot_name)
msg.attach(part)


server = smtplib.SMTP("smtp.gmail.com:587")
server.ehlo()
server.starttls()
server.login("assemblyshine@gmail.com", "rshiny@2020")

server.sendmail(msg['From'], emaillist , msg.as_string())
