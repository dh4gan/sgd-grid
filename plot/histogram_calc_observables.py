# Written by Duncan Forgan, 17/1/2013
# Code reads output from FORTRAN code selfgravdisc_modelgrid
# Creates a histogram of the requested variables

import numpy as np
import matplotlib.pyplot as plt
import filefinder as ff

# Set up tuples and dictionaries

variablekeys = ("q", "mtot", "flux_nu","fluxtot", "mflux", "mtot_mflux")
variablenames = ("Disc to Star Mass Ratio","Actual Disc Mass",r"$F_{\nu}$",r"$F_{\rm tot}$","Estimated Disc Mass", "Estimated Mass / True Mass")
variablecolumns = (2,3,4,5,6,7)


namedict = {}
coldict = {}

for i in range(len(variablekeys)):
    namedict[variablekeys[i]] = variablenames[i]
    coldict[variablekeys[i]] = variablecolumns[i]


# Read in filename
inputfile = ff.find_local_input_files('*.observe')

# Decide which variable to plot contours for

print "Which variables are to be binned?"

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
    
nbins = input("How many bins?   ")
cumulative = raw_input("Cumulative distribution? (y/n)  ")
    
print "Reading File ", inputfile

# Read header first
f = open(inputfile, 'r')

header = f.readline()
headernums = header.split()

nrad = int(headernums[0])
nmdot = int(headernums[1])
nu = float(headernums[2])
lam = float(headernums[3])
mstar = float(headernums[8])

print 'Observables data at wavelength ',lam, ' microns'

if(nu > 1.0e9):
    print 'Frequency: ',nu/1.0e9, ' GHz'
elif(nu > 1.0e6):
    print 'Frequency: ',nu/1.0e9, ' MHz'

# Now read rest of file

data = np.genfromtxt(inputfile, skiprows=1)
print "File Read"

# Now loop over choices

for i in range(len(choices)):
    # Extract data column and reshape
    print "Plotting ",namedict[choices[i]]
    plotdata = data[:,columns[i]]
    plotdata = plotdata[:][plotdata!=0.0] # Delete all data=zero
    plotdata = plotdata[:][np.abs(plotdata)<1.0e25] # Delete all nonsensical data!
    
    # If there is no data to plot, skip it    
    if not plotdata.size:
        print "Skipping plot: no data"        
        continue
    
    # Make plot
    
    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    ax.set_ylabel("Relative Frequency")
    ax.set_xlabel(namedict[choices[i]])
    
    if cumulative=="y":
        hist = ax.hist(plotdata,bins=nbins, normed=True,cumulative = True )
    else:
        hist = ax.hist(plotdata,bins=nbins, normed=True,cumulative = False)
    plt.title(namedict[choices[i]]+r', $\lambda=$'+str(lam) + r'$\mu m$')
    outputfile = "hist"+choices[i]+'.ps'
    plt.savefig(outputfile, format = 'ps')
