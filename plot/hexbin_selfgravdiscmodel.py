# Written by Duncan Forgan, 17/1/2013
# Code reads output from FORTRAN code selfgravdisc_modelgrid
# Creates a hexbin plot of the requested variable in mdot-r space

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys

import filefinder as ff

# Set up tuples and dictionaries

variablekeys = ("q","sigma", "cs","omega","temp", "betac","alpha","mjeans","ljeans","rhill","h","tirr", "csirr", "Qirr", "gamma_J", "gamma_Q")
variablenames = ("Disc to Star Mass Ratio","Surface Density","Sound Speed","Angular Velocity","Temperature","Beta_c", "Alpha",r"$M_{\rm Jeans}$ ($M_{\rm Jup}$)","Jeans length", "Hill Radius", "H/R", "Tirr", "cs_irr","Qirr", "Gamma_J", "Gamma_Q")
variablecolumns = (2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)


namedict = {}
coldict = {}

for i in range(len(variablekeys)):
    namedict[variablekeys[i]] = variablenames[i]
    coldict[variablekeys[i]] = variablecolumns[i]


# Read in filename
inputfile = ff.find_local_input_files('*.sgdmodel')
contourchoice = raw_input("Add contour lines? (y/n) ")

add_contour=False
if("y" in contourchoice or "Y" in contourchoice):
    add_contour=True

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
    
print "Reading File ", inputfile

# Read header first
f = open(inputfile, 'r')

header = f.readline()
headernums = header.split()

nrad = int(headernums[0])
nmdot = int(headernums[1])
mstar = float(headernums[6])

# Now read rest of file

data = np.genfromtxt(inputfile, skiprows=1)

print "File Read"

# Pull radius and mdot data, and reshape
rad = data[:,0]
mdot = data[:,1]

# Now loop over choices

for i in range(len(choices)):
    # Extract data column and reshape
    print "Plotting ",namedict[choices[i]]
    
    plotdata = data[:,columns[i]]
    contourrad = rad.reshape(nmdot,nrad)
    contourmdot = mdot.reshape(nmdot,nrad)
    contourdata = plotdata.reshape(nmdot,nrad)
    indices = plotdata[:]>1e-40

    print contourdata.shape
    # Delete junk data

    #indices = plotdata[:]>1.0e-50
    
    radplot = rad[indices]
    mdotplot = mdot[indices]
    plotdata = plotdata[indices] # Delete all nonsensical data!

    
    #if choices[i]=='mjeans':
    #    plotdata = plotdata*0.000954        
        
    if choices[i]=='h':
        plotdata = plotdata/radplot
    
    plotmin = np.amin(plotdata)
    plotmax = np.amax(plotdata)
    
    
    
    print 'min, max: ',plotmin,plotmax
    #if choices[i]=='mjeans':
        #plotmin = 0.08
        #plotmax = 0.35
    
    #if choices[i]=='q':
    #    plotmin = 0.1
    #    plotmax = 1.0 
        
    if choices[i]=='temp':
        plotmin = 50.0           
    
    print plotmin, plotmax
    
    if(plotmin == plotmax):
        print "Skipping plot: plotmin=plotmax"
        print plotmin, plotmax
        continue
    
    
    # Make plot
    
    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    ax.set_ylabel(r"Accretion Rate, $\mathrm{(M_{\odot} yr^{-1})}$", fontsize=20)
    ax.set_xlabel(r"$r_{out}$ (AU)", fontsize = 20)
    ax.set_yscale('log')
    plt.hexbin(radplot,mdotplot,C=plotdata,gridsize = int(nrad*0.25), vmin = plotmin, vmax = plotmax, yscale='log',mincnt = 1,cmap='Blues')


   
    #if(choices[i]=='mjeans'):
	#ax.set_ylim(3e-5,1e-3)    

    # Add a hatched background where model does not return a value
    ax.set_axis_bgcolor('gray')
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xy = (xmin,ymin)
    width = xmax - xmin
    height = ymax -ymin
    p = patches.Rectangle(xy, width, height, hatch='x', fill=None, zorder=-10)
    ax.add_patch(p)


    cb = plt.colorbar()
    cb.set_label(namedict[choices[i]], fontsize=20)
    
    # Add contours
    if(add_contour):
         cs = plt.contour(contourrad,contourmdot,contourdata,colors='white')
         plt.clabel(cs,cs.levels,inline=True)


    #plt.title(namedict[choices[i]])
    outputfile = "mdot_r_"+choices[i]+'.png'
    plt.savefig(outputfile, format = 'png')
