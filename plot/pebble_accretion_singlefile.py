# Written 7/1/17 by dh4gan
# Plots output data from calc_pebble_accretion
# Reads individual files for a specific gas accretion rate

import numpy as np
import matplotlib.pyplot as plt
import filefinder as ff


# Set up tuples and dictionaries

variablekeys = ("rpeb","rdotpeb","mdotpeb","rmin_stream","rmax_stream","mcross", "planetmdotpeb","eff")
variablenames = ("$r_{\rm peb}$", "$\dot{r}_{\rm peb}$","$\dot{M}_{\rm peb}\,(M_{\rm } \, \rm{yr}^{-1}$)","$r_{\rm min, stream}$","$r_{\rm max,stream}$","$M_{\rm cross} (M_{\rm Jup})$","\dot{M}_{pl}$", "\epsilon")
variablecolumns = range(8)

namedict = {}
coldict = {}

for i in range(len(variablekeys)):
    namedict[variablekeys[i]] = variablenames[i]
    coldict[variablekeys[i]] = variablecolumns[i]

# Find pebble accretion file

filename = ff.find_local_input_files('*.mdotpebble.*')

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
mdotgas = float(headernums[4])
mstar = float(headernums[5])
tstop = float(headernums[6])

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
    ax.set_xlabel(namedict["rpeb"])
    ax.set_yscale('log')
    ax.plot(data[:,coldict["rpeb"]), data[:,coldict[choices[i]])

    outputfile = filename+choices[i]+".png"

    fig1.savefig(outputfile)


# Plot inner and outer streaming instability limits as well

print "Also plotting inner and outer streaming instability boundaries for this accretion rate"

fig1 = plt.figure()

ax1 = fig1.add_subplot(111)

ax1.plot(data[:,coldict["rpeb"]],data[:,coldict["rmin"])
ax1.plot(data[:,coldict["rpeb"]],data[:,coldict["rmax"])
ax1.set_xlabel(variablenames["rpeb"])
ax1.set_ylabel("Streaming Instability Boundaries")

outputfile = filename+"stream_limits.png"
fig1.savefig(filename)

