#!/usr/bin/env python3.6

import math
import matplotlib
#matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import cm
from seppyio import *
from matplotlib.patches import FancyBboxPatch 
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import LinearLocator

filename=searchArgv('input') 

n1,o1,d1=from_header('input','n1','o1','d1')

n1=int(n1)

dat=read('input',(n1))

# compute data stat
datRMS = np.linalg.norm(dat)/math.sqrt(dat.size)
datMIN = np.min(dat)
datMAX = np.max(dat)

# read parameters
zmin=searchArgv('xmin')
if zmin!=False:
    zmin=float(zmin)
else:
    zmin=o1

zmax=searchArgv('xmax')
if zmax!=False:
    zmax=float(zmax)
else:
    zmax=o1+d1*n1

izmin=int((zmin-o1)/d1)
izmax=int((zmax-o1)/d1)

izmin=max(0,izmin)
izmax=min(n1,izmax)

dat=dat[izmin:izmax+1]
#dat=np.flip(dat,axis=0)

xlabel=searchArgv('xlabel')
if xlabel!=False:
    xlabel=str(xlabel)
else:
    xlabel='X (m)'

zlabel=searchArgv('zlabel')
if zlabel!=False:
    zlabel=str(zlabel)
else:
    zlabel='Time (s)'

gain=searchArgv('gain')
if gain!=False:
    gain=float(gain)
else:
    gain=1
dat=gain*dat


min=searchArgv('min')
if min!=False:
    min=float(min)
else:
    min=-datRMS

max=searchArgv('max')
if max!=False:
    max=float(max)
else:
    max=datRMS


normalize=searchArgv('normalize')
if normalize!=False:
    normalize=str(normalize)

if normalize=='minmax':
    min=datMIN
    max=datMAX
elif normalize=='rms':
    min=-datRMS
    max=datRMS

zreverse=searchArgv('zreverse')
if zreverse==False:
    zreverse=0
elif zreverse=='0':
    zreverse=0
else:
    zreverse=1

format=searchArgv('format')
if format!=False:
    format=str(format)
else:
    format='png'

title=searchArgv('title')
if title==False:
    title=''

fontsize=searchArgv('fontsize')
if fontsize!=False:
    fontsize=int(fontsize)

zlogscale=searchArgv('zlogscale')
if zlogscale==False:
    zlogscale=0
elif zlogscale=='0':
    zlogscale=0
else:
    zlogscale=1

xsize=searchArgv('xsize')
if xsize==False:
    xsize=6
else:
    xsize=float(xsize)

zsize=searchArgv('zsize')
if zsize==False:
    zsize=4
else:
    zsize=float(zsize)

transpose=searchArgv('transpose')
if transpose==False:
    transpose=0
elif transpose=='0':
    transpose=0
else:
    transpose=1


x=np.arange(zmin, zmax, d1)

fig, ax = plt.subplots(figsize=(xsize, zsize))
if transpose==0:
    ax.plot(x,dat)
    plt.xlabel(xlabel)
    plt.ylabel(zlabel)
    plt.xlim(zmin, zmax)
    plt.ylim(min,max)
else:
    ax.plot(dat,x)
    plt.ylabel(xlabel)
    plt.xlabel(zlabel)
    plt.ylim(zmin, zmax)
    plt.xlim(min,max)

plt.title(title)
#plt.xticks(np.arange(zmin,zmax,d1))


if zreverse==1:
    plt.gca().invert_yaxis()

if zlogscale==1:
    plt.yscale('log')

out=searchArgv('output')
if out==False: 
    plt.show()
else:
    plt.savefig(out,bbox_inches='tight',format=format)