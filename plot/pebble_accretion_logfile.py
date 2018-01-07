# Written 7/1/17 by dh4gan
# Plots output data from calc_pebble_accretion
# Reads the log file, computed over all gas accretion rates

import numpy as np
import matplotlib.pyplot as plt
import filefinder as ff


# Set up tuples and dictionaries

variablekeys = ("mdotgas","rpebmax","tpebmax","planetmdotpebmax","mcrossmax", "effmax")
variablenames = ("$\dot{M}_{\rm gas}$","$r_{\rm peb,max}$", "$t_{\rm peb,max}$","$\dot{M}_{pl,max}$","$M_{\rm cross,max} (M_{\rm Jup})$","\epsilon_{\rm max}")
variablecolumns = range(6)

namedict = {}
coldict = {}

for i in range(len(variablekeys)):
    namedict[variablekeys[i]] = variablenames[i]
    coldict[variablekeys[i]] = variablecolumns[i]

# Find pebble accretion file

filename = ff.find_local_input_files('*.mdotpebblelog')

# Decide which variable to plot

print "Which variables are to be plotted?"

for i in range(len(variablekeys)):
    print variablekeys[i],":\t \t \t", namedict[variablekeys[i]]

print "all : \t \t \t Plots all variables"
keyword = raw_input("Enter appropriate keywords separated by spaces:   ")


# If all selected, then generate all keywords automatically

if "all" in keyword:
    choices = variablekeys
else:
    # Otherwise, parse keyword string into individual choices
    choices = keyword.split()

columns = []

# Determine variable columns from keywords
for word in choices:
    columns.append(coldict[word])

# Open file and read data

print "The following columns are to be plotted:"
for i in range(len(choices)):
    print namedict[choices[i]],columns[i]
    
print "Reading File ", inputfile

# Read header first
f = open(inputfile, 'r')

header = f.readline()
headernums = header.split()

nrad = int(headernums[0])
nmdot = int(headernums[1])
mstar = float(headernums[6])
tstop = float(headernums[7])

# Now read rest of file

data = np.genfromtxt(inputfile, skiprows=1)

print "File Read"


# Now loop over choices

for i in range(len(choices)):
    # Extract data column and reshape
    print "Plotting ",namedict[choices[i]]

    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    ax.set_ylabel(namedict[choices[i])
    ax.set_xlabel(namedict["mdotgas"])
    ax.set_yscale('log')
    ax.plot(data[:,coldict["mdotgas"]), data[:,coldict[choices[i]])

    outputfile = filename+choices[i]+".png"

    fig1.savefig(outputfile)

