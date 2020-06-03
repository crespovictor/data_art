# -------------------------------------------- # 
# ----------------- setup -------------------- #
# -------------------------------------------- # 

# import required libraries

from Bio import SeqIO
import pandas as pd
import holoviews as hv
from holoviews import opts, dim

# prepare holoviews canvas

hv.extension('bokeh')
hv.output(size=300)

# auxiliar variables used to process data and create relationships

lista = []
dictionary = {"AA":0,"AC":0,"AG":0,"AT":0,  # all possible relationships of nucletids
              "CA":0,"CC":0,"CG":0,"CT":0,
              "GA":0,"GC":0,"GG":0,"GT":0,
              "TA":0,"TC":0,"TG":0,"TT":0}
dict_aux = {"AA":0,"AC":1,"AG":2,"AT":3,   # indexes of relationships used to add weights to
            "CA":4,"CC":5,"CG":6,"CT":7,   # data
            "GA":8,"GC":9,"GG":10,"GT":11,
            "TA":12,"TC":13,"TG":14,"TT":15}


# storage of links and weights used to create the diagram

data = [["A", "A", 0],["A", "C", 0],["A", "G", 0],["A", "T", 0],
        ["C", "A", 0],["C", "C", 0],["C", "G", 0],["C", "T", 0],
        ["G", "A", 0],["G", "C", 0],["G", "G", 0],["G", "T", 0],
        ["T", "A", 0],["T", "C", 0],["T", "G", 0],["T", "T", 0]]



# store the name of the FASTA file in a variable

fasta_file = "EPI_ISL_419959.fasta"

# open the FASTA file
# this will store all the data from the file
# the outout will be something like this:
#
# ID: hCoV-19/Australia/VIC266/2020|EPI_ISL_419959
# Name: hCoV-19/Australia/VIC266/2020|EPI_ISL_419959
# Description: hCoV-19/Australia/VIC266/2020|EPI_ISL_419959
# Number of features: 0
# Seq('ACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGC...AAA', SingleLetterAlphabet())
#
# we are only interested in the seq part of the file

mito = SeqIO.read(fasta_file,"fasta")


# store the sequence in a string

sequence = mito.seq

# -------------------------------------------- # 
# ------------ data wranglings --------------- #
# -------------------------------------------- # 

# transform the string into a list of characters
for letter in sequence:
    lista.append(letter)

# store the number of relationships occurred (i.e, AA, AC, AT, etc)    
for i in range(len(lista)-1):
    dictionary[str(lista[i]+lista[i+1])] = dictionary[str(lista[i]+lista[i+1])] +1
    
# create the dataset
for i in range(len(data)):
    data[i][2] = dictionary[data[i][0]+data[i][1]]

# -------------------------------------------- # 
# ------------ create diagram ---------------- #
# -------------------------------------------- # 

# plot the chart using holoviews
hm = hv.HeatMap(data).sort()

# this line creates a diagram with labels, comment it if you don't want labels
# this is useful if you want to understand what the colors mean
hm.opts(colorbar = True, xlabel = "Source", ylabel = "Target")


# this line creates a diagram without labels or anotations
# this is useful if you want to create a more artsy visualizations
#hm.opts(xaxis = None, yaxis = None)