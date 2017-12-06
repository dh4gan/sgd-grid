# Written by Duncan Forgan, 17/1/2013
# Code reads output from FORTRAN code selfgravdisc_modelgrid --> calc_observables
# Creates a 2D contour plot of the requested variables in mdot-r space

import numpy as np
import matplotlib.pyplot as plt
import filefinder as ff

# Set up tuples and dictionaries

variablekeys = ("q", "mtot", "flux_nu","fluxtot", "mflux", "mtot_mflux")
variablenames = ("Disc to Star Mass Ratio","Actual Disc Mass",r"$F{r)$",r"Total Disc Flux","Estimated Disc Mass", "Estimated Mass / True Mass")
variablecolumns = (2,3,4,5,6,7)

xobs = 100.0
yobs = 6.6e-7

namedict = {}
coldict = {}

for i in range(len(variablekeys)):
    namedict[variablekeys[i]] = variablenames[i]
    coldict[variablekeys[i]] = variablecolumns[i]


# Read in filename
inputfile = ff.find_local_input_files('*.observe')

# Decide which variable to plot contours for

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
    
ncontour = input("How many contour lines? "   )
    
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

lamunit = r'$\mu $m'

# Change wavelength unit to mm if necessary
if(lam < 1.0):
    lamunit = r'nm'
    lam = lam*1.0e3
    print 'Equivalently, wavelength '+str(lam) +'  '+lamunit

if(lam > 1.0e3):
    lamunit = r'mm'
    lam = lam/1.0e3
    print 'Equivalently, wavelength '+str(lam) +'  '+lamunit

    if(lam > 1.0e2):
        lamunit = r'cm'
        lam = lam/1.0e2
        print 'Equivalently, wavelength '+str(lam) +'  '+lamunit

if(nu > 1.0e9):
    print 'Frequency: ',nu/1.0e9, ' GHz'
elif(nu > 1.0e6):
    print 'Frequency: ',nu/1.0e9, ' MHz'

# Now read rest of file

data = np.genfromtxt(inputfile, skiprows=1)

print "File Read"
# Pull radius and mdot data, and reshape
rad = data[:,0]
mdot = data[:,1]

rad = rad.reshape(nrad,nmdot)
mdot = mdot.reshape(nrad,nmdot)

# Now loop over choices

for i in range(len(choices)):
    # Extract data column and reshape
    print "Plotting ",namedict[choices[i]]
    
    plotdata = data[:,columns[i]]    
    plotdata = plotdata[:][plotdata<1.0e25] # Delete all nonsensical data!
    plotdata = plotdata.reshape(nrad,nmdot)    
    plotmin = np.amin(plotdata)
    plotmax = np.amax(plotdata)
    
    if(plotmin == plotmax):
        print "Skipping plot: plotmin=plotmax"
        print plotmin, plotmax
        continue
    
    
    # Make plot
    
    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    ax.set_ylabel(r"Accretion Rate, $\mathrm{(M_{\odot} yr^{-1})}$")
    ax.set_xlabel(r"r (AU)")
    ax.set_yscale('log')
    cont = ax.contour(rad,mdot,plotdata,ncontour,colors='black')
    obs = ax.scatter(xobs,yobs,marker='x',s=100,c='red')
    if choices[i] !='mtot':
        plt.title(namedict[choices[i]]+r', $\lambda=$'+str(lam) + ' '+lamunit)
    else:
        plt.title(namedict[choices[i]])
    
    #plt.clabel(cont,cont.levels[1::2], inline=1, fontsize=10)
    
    manual_locations = []
    dr = 150.0/ncontour
    for c in range(ncontour):        
        manual_locations.append([c*dr,1e-7])    
            
    plt.clabel(cont,cont.levels[1::2,], manual=manual_locations)

    outputfile = "mdot_r_"+choices[i]+'.ps'
    plt.savefig(outputfile, format = 'ps')
