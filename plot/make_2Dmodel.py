# Written by Duncan Forgan, 4/10/2016
# Code reads output from FORTRAN code selfgravdisc_modelgrid (1D, radial model)
# Creates a 2D (r,phi) model

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import filefinder as ff

def find_nearest(array,value):
    return np.abs(array-value).argmin()

# Set up tuples and dictionaries

variablekeys = ("q","sigma", "cs","omega","temp", "betac","alpha","mjeans","ljeans","rhill","h","tirr", "csirr", "Qirr", "gamma_J", "gamma_Q")
variablenames = ("Disc to Star Mass Ratio","Surface Density","Sound Speed","Angular Velocity","Temperature","Beta_c", "Alpha","$M_{Jeans}/M_{\odot}$","Jeans length", "Hill Radius", "H/R", "Tirr", "cs_irr","Qirr", "Gamma_J", "Gamma_Q")
variablecolumns = (2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)

ncol=17
namedict = {}
coldict = {}

udist = 1.496e8 # AU --> km

for i in range(len(variablekeys)):
    namedict[variablekeys[i]] = variablenames[i]
    coldict[variablekeys[i]] = variablecolumns[i]


# Now define logarithmic spiral parameters (see Hall et al 2016, MNRAS)

a_spiral = 13.5
b_spiral = 0.38
m_spiral = 2

# transport locality

#eta_trans = 0.75
eta_trans = 0.0



nx = 200
ny = 200

# Read in filename
inputfile = ff.find_local_input_files('*.model')


    
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

# Pull radius, mdot and sigma columns
raddat = data[:,0]
mdotdat = data[:,1]
sigma = data[:,coldict['sigma']]
cs = data[:,coldict['cs']]

# Find the unique values in the array
mdot = np.unique(mdotdat)
rad = np.unique(raddat)

rmin = rad[0]
rmax = rad[-1]

# Which model is going to be made into 2D?
#mdotchoice = input("What accretion rate to use? \n min="+str(mdot[0])+", max="+str(mdot[-1])+ "\n")
#rchoice = input("What is the outer radius? \n min="+str(rmin)+", max="+str(rmax)+ "\n")

mdotchoice = 1.0e-4
rchoice = 1200.0

# Find model corresponding to these parameters

imdot = find_nearest(mdot,mdotchoice)
irad = find_nearest(rad,rchoice)

print "Located disc model at mdot, radius indices: ", imdot, irad

mdotchoice = mdot[imdot]
rchoice = rad[irad]

print "mdot: ", mdotchoice
print "rad: ", rchoice

# Extract sigma
# Find first instance of mdot in the column

firstindex = np.where(mdotdat==mdotchoice)[0][0]
sigma1D = data[firstindex:firstindex+nrad,coldict['sigma']]
alpha1D = data[firstindex:firstindex+nrad,coldict['alpha']]
omega1D = data[firstindex:firstindex+nrad,coldict['omega']]
cs1D = data[firstindex:firstindex+nrad,coldict['cs']]



fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(rad, sigma1D, label= 'sigma')
ax1.plot(rad, alpha1D, label= 'alpha')
ax1.plot(rad, omega1D, label= 'omega')
ax1.plot(rad,cs1D,label='cs')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.legend()


# Define 2D arrays to hold sigma, velocity

x = np.linspace(-rmax,rmax,num=nx)
y = np.linspace(-rmax,rmax,num=ny)

sigma2D = np.zeros((nx,ny))
cs2D = np.zeros((nx,ny))
vel2D = np.zeros((nx,ny))
vx2D = np.zeros((nx,ny))
vy2D = np.zeros((nx,ny))

# Now construct 2D array

for ix in range(nx):
    for iy in range(ny):
        
        r = np.sqrt(x[ix]*x[ix]+ y[iy]*y[iy])
        phi = np.arctan2(y[iy],x[ix])
        
        
        
        phi_spiral = np.mod(np.log(r/a_spiral)/b_spiral, 2.0*np.pi)
        
        ir = find_nearest(rad,r)
        
        cs2D[ix,iy] = cs1D[ir]
        
        #dsigma = np.sqrt(alpha1D[ir])
        dsigma = 0.0
        dsigma = -dsigma*np.cos(m_spiral*(phi_spiral-phi))
                
        dvel = -(eta_trans)*np.cos(m_spiral*(phi_spiral-phi))
        
        sigma2D[ix,iy] = sigma1D[ir]*(1.0+dsigma)
        vel2D[ix,iy] = omega1D[ir]*r*udist*(1.0+dvel)
        
        vx2D[ix,iy] = vel2D[ix,iy] * np.sin(phi)
        vy2D[ix,iy] = vel2D[ix,iy] * np.cos(phi)
        
        
# Do a 2D plot?

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
plt.title(r'Surface Density $(g cm^-2)$')
plt.pcolor(x,y,sigma2D, norm=LogNorm(vmin=sigma2D.min(), vmax=sigma2D.max()), cmap = 'PuBu_r')
plt.colorbar()


fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
plt.title(r'Velocity ')
plt.pcolor(x,y,vel2D, norm=LogNorm(vmin=vel2D.min(), vmax=vel2D.max()), cmap = 'PuBu')
plt.colorbar()
plt.show()



print 'Calculations complete, writing to file'
# Write the 2D data to file

outputfile = inputfile+'.2D'

f_obj = open(outputfile, 'w')
line = str(mstar)+ "\t" + str(mdotchoice) +"\t"+ str(rmin) +"\t"+ str(rmax) + "\n"
f_obj.write(line)

for ix in range(nx):
    for iy in range(ny):
        line = str(x[ix])+ "\t"+str(y[iy])+"\t"+str(sigma2D[ix,iy])+"\t"+str(cs2D[ix,iy])+"\t"+str(vx2D[ix,iy])+"\t"+str(vy2D[ix,iy])+"\n"
        f_obj.write(line)

f_obj.close()
