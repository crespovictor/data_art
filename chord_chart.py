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

data = [{'source': 1, 'target': 1, 'value': 0},{'source': 1, 'target': 2, 'value': 0},
        {'source': 1, 'target': 3, 'value': 0},{'source': 1, 'target': 4, 'value': 0},
        {'source': 2, 'target': 1, 'value': 0},{'source': 2, 'target': 2, 'value': 0},
        {'source': 2, 'target': 3, 'value': 0},{'source': 2, 'target': 4, 'value': 0},
        {'source': 3, 'target': 1, 'value': 0},{'source': 3, 'target': 2, 'value': 0},
        {'source': 3, 'target': 3, 'value': 0},{'source': 3, 'target': 4, 'value': 0},
        {'source': 4, 'target': 1, 'value': 0},{'source': 4, 'target': 2, 'value': 0},
        {'source': 4, 'target': 3, 'value': 0},{'source': 4, 'target': 4, 'value': 0}]



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

# add values to the datased used for the actual diagram     
for k in dict_aux:
    data[dict_aux[k]]["value"] = dictionary[k]
    

# -------------------------------------------- # 
# ------------ create diagram ---------------- #
# -------------------------------------------- # 

# transform data to a Pandas dataset
links = pd.DataFrame(data)
nodos = [{"name":"Adenine", "index":1},{"name":"Cytosine", "index":2},
         {"name":"Guanine", "index":3},{"name":"Thymine", "index":4}]
nodes = hv.Dataset(pd.DataFrame(nodos), 'index')

# generate the diagram with given options
chord = hv.Chord((links, nodes))
chord.opts(
    opts.Chord(cmap='Category20', edge_cmap='Category20', edge_color=dim('source').str(), 
               labels="name", node_color=dim('index').str()))
