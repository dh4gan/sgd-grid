# Written 7/1/17 by dh4gan
# Plots output data from calc_pebble_accretion
# Reads individual files for a specific gas accretion rate

import numpy as np
import matplotlib.pyplot as plt
import filefinder as ff

plt.rc('axes',labelsize=18)

# Set up tuples and dictionaries

variablekeys = ("rpeb","tpeb","rdotpeb","mdotpeb","rmin_stream","rmax_stream","mcross","mjeans", "planetmdotpeb","planeteff")
variablenames = (r"$r_{\rm peb}$ (AU)", r"$t_{\rm peb}$ (yr)",r"$\dot{r}_{\rm peb}$ (AU yr$^{-1}$)",r"$\dot{M}_{\rm peb}\,(M_{\rm } \, \rm{yr}^{-1}$)",r"$r_{\rm min, stream}$",r"$r_{\rm max,stream}$",r"$M_{\rm cross} (M_{\rm Jup})$",r"$M_{\rm jeans} (M_{\rm Jup})$ ",r"\dot{M}_{pl}$", "\epsilon")
variablecolumns = range(len(variablekeys))

namedict = {}
coldict = {}

for i in range(len(variablekeys)):
    namedict[variablekeys[i]] = variablenames[i]
    coldict[variablekeys[i]] = variablecolumns[i]

# Find pebble accretion file

filename = ff.find_sorted_local_input_files('*.pebble')

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
    
print "Reading File ", filename

# Read header first
f = open(filename, 'r')

header = f.readline()
headernums = header.split()

nrad = int(headernums[0])
nmdot = int(headernums[1])
mdotgas = float(headernums[4])
mstar = float(headernums[5])
tstop = float(headernums[6])

# Now read rest of file

data = np.genfromtxt(filename, skiprows=1)

print "File Read"


# Now loop over choices

for i in range(len(choices)):
    # Extract data column and reshape
    print "Plotting ",namedict[choices[i]]

    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    ax.set_ylabel(namedict[choices[i]])
    ax.set_xlabel(namedict["rpeb"])
    ax.set_yscale('log')
    ax.plot(data[:,coldict["rpeb"]], data[:,coldict[choices[i]]])

    outputfile = choices[i]+"vs_rpeb"+filename+".png"

    fig1.savefig(outputfile)


# Plot inner and outer streaming instability limits as well

print "Also plotting inner and outer streaming instability boundaries for this accretion rate"

fig1 = plt.figure()

ax1 = fig1.add_subplot(111)

ax1.plot(data[:,coldict["rpeb"]],data[:,coldict["rmin_stream"]], color='blue')
ax1.plot(data[:,coldict["rpeb"]],data[:,coldict["rmax_stream"]],color = 'green')
ax1.set_xlabel(namedict["rpeb"])
ax1.set_ylabel("Streaming Instability Boundaries (AU)")

outputfile = "stream_limits"+filename+".png"
fig1.savefig(outputfile)

