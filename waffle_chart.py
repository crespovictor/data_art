# -------------------------------------------- # 
# ----------------- setup -------------------- #
# -------------------------------------------- # 

# import required libraries
from Bio import SeqIO
import math
import matplotlib.pyplot as plt
from pywaffle import Waffle

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

# get percentages of each nucleotid
# TODO: develop a function to get the right amounts comparing them and avoid hardcoding
a = math.ceil(sequence.count("A")/len(sequence) * 100)-1 # harcoded to substract one.
c = math.ceil(sequence.count("C")/len(sequence) * 100)
g = math.ceil(sequence.count("G")/len(sequence) * 100)
t = math.ceil(sequence.count("T")/len(sequence) * 100)

# with the percentages, create the dataset to use
data = {"A": a, "C": c, "G": g, "T":t}

# -------------------------------------------- # 
# ------------ create diagram ---------------- #
# -------------------------------------------- # 

# plot the chart using matplotlib
fig = plt.figure(
    FigureClass=Waffle, 
    rows=10, 
    columns=10,
    values=data, 
    legend={'loc': 'upper left', 'bbox_to_anchor': (1, 1)},
    figsize=(10, 10)  # figsize is a parameter of matplotlib.pyplot.figure
)
plt.show()