# Written by Duncan Forgan, 7/1/2017
# Code reads output from calc_pebble_accretion
# Creates a hexbin plot of the requested variable in mdotgas-rpeb space

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys

#sys.path.append('/data/dhf3/programs/python/filefinder')
import filefinder as ff

# Set up tuples and dictionaries


variablekeys = ("grainsize", "tstop","tstopratio","maxgrainsize","tpeb","rdotpeb","mdotpeb","rmin_stream","rmax_stream","mcross","mjeans", "planetmdotpeb","planeteff")
variablenames = (r"$s$ (cm)", r"$\tau_s$",r"$\tau_s/\tau_{s,max}$",r"$s_{\rm max}$ cm",r"$t_{\rm peb}$ (yr)",r"$\dot{r}_{\rm peb}$ (AU yr$^{-1}$)",r"$\dot{M}_{\rm peb}\,(M_{\rm } \, \rm{yr}^{-1}$)",r"$r_{\rm min, stream}$",r"$r_{\rm max,stream}$",r"$M_{\rm cross} (M_{\rm Jup})$",r"$M_{\rm jeans} (M_{\rm Jup})$ ",r"$\dot{M}_{pl}$", r"$\epsilon$")
variablecolumns = range(2,len(variablekeys)+2)

namedict = {}
coldict = {}

print variablecolumns

for i in variablecolumns:

    namedict[variablekeys[i-2]] = variablenames[i-2]
    coldict[variablekeys[i-2]] = i
    print i, variablekeys[i-2], variablenames[i-2]


# Read in filename
inputfile = ff.find_sorted_local_input_files('*.pebble')

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
#mstar = float(headernums[6])

# Now read rest of file

data = np.genfromtxt(inputfile, skiprows=1)

print "File Read"

# Pull radius and mdot data, and reshape
rad = data[:,1]
mdot = data[:,0]

# Now loop over choices

for i in range(len(choices)):
    # Extract data column and reshape
    print "Plotting ",namedict[choices[i]]
    
    plotdata = data[:,columns[i]]
    #indices = plotdata[:]>1e-40
    
    # Delete junk data

    #indices = plotdata[:]>1.0e-50
    radplot = rad
    mdotplot = mdot
    #radplot = rad[indices]
    #mdotplot = mdot[indices]
    #plotdata = plotdata[indices] # Delete all nonsensical data!
    
    plotmin = np.amin(plotdata)
    plotmax = np.amax(plotdata)
    
    
    print 'min, max: ',plotmin,plotmax

    if choices[i]=='tstopratio':
        plotmin = 0.0
        plotmax = 1.0
#    if choices[i]=='mjeans':
        #plotmin = 0.08
#        plotmax = 0.35
    
#    if choices[i]=='q':
#        plotmin = 0.1
#        plotmax = 1.0 
        
#    if choices[i]=='temp':
#        plotmin = 50.0           
    
#    print plotmin, plotmax
    
    if(plotmin == plotmax):
        print "Skipping plot: plotmin=plotmax"
        print plotmin, plotmax
        continue
    
    
    # Make plot
    
    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    ax.set_ylabel(r"Gas Accretion Rate, $\mathrm{(M_{\odot} yr^{-1})}$", fontsize=20)
    ax.set_xlabel(r"$r_{out}$ (AU)", fontsize = 20)
    ax.set_yscale('log')
    plt.hexbin(radplot,mdotplot,C=plotdata,gridsize = int(nrad*0.25), vmin = plotmin, vmax = plotmax, yscale='log',mincnt = 1,cmap='Blues')

    if(choices[i]=='mjeans'):
	ax.set_ylim(3e-5,1e-3)    

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
    #plt.title(namedict[choices[i]])
    outputfile = "mdot_r_"+choices[i]+'.png'
    plt.savefig(outputfile, format = 'png')
