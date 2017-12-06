# Written by Duncan Forgan, 17/1/2013
# Code reads output from FORTRAN code selfgravdisc_modelgrid --> calc_observables
# Creates a 2D contour plot of the requested variables in mdot-r space

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys
sys.path.append('/disk1/dhf/programs/python/filefinder')
import filefinder.localfiles as ff

# Set up tuples and dictionaries

variablekeys = ("q", "mtot", "flux_nu","fluxtot", "mflux", "mtot_mflux")
variablenames = (r"Disc to Star Mass Ratio",r"Actual Disc Mass ($M_{\odot}$)",r"$F{r)$",r"Total Disc Flux (mJy)","Estimated Disc Mass ($M_{\odot}$) ", "Estimated Mass / True Mass")
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

# Now loop over choices

for i in range(len(choices)):
    # Extract data column and reshape
    print "Plotting ",namedict[choices[i]]
    
    plotdata = data[:,columns[i]]           
    indices = plotdata[:]>1.0e-40
    
    # Delete junk data    
    
    radplot = rad[indices]
    mdotplot = mdot[indices]
    plotdata = plotdata[indices] # Delete all nonsensical data!             
    
    plotmin = np.amin(plotdata)
    plotmax = np.amax(plotdata)
    
    if(plotmin == plotmax):
        print "Skipping plot: plotmin=plotmax"
        print plotmin, plotmax
        continue
                
    # Make plot
    
    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    ax.set_ylabel(r"Accretion Rate, $\mathrm{(M_{\odot} yr^{-1})}$",fontsize=20)
    ax.set_xlabel(r"$r_{out}$ (AU)",fontsize=20)
    ax.set_yscale('log')
    ax.set_axis_bgcolor('gray')
    #ax.set_ylim(np.amin(mdot), 1.0e-4)
    plt.hexbin(radplot,mdotplot,C=plotdata,gridsize = int(nrad*0.25), vmin = plotmin, vmax = plotmax, yscale='log',mincnt = 1,cmap='Blues')

    # Add a hatched background where model does not return a value
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xy = (xmin,ymin)
    width = xmax - xmin
    height = ymax -ymin
    p = patches.Rectangle(xy, width, height, hatch='x', fill=None, zorder=-10)
    ax.add_patch(p)
    cb = plt.colorbar()
        
    if choices[i] !='mtot':
        cb.set_label(namedict[choices[i]]+r', $\lambda=$'+str(lam) + ' '+lamunit, fontsize=20)
    else:
        cb.set_label(namedict[choices[i]], fontsize=20)

    outputfile = "mdot_r_"+choices[i]+'.png'
    plt.savefig(outputfile, format = 'png')
    
