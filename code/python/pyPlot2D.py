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

filename=searchArgv('input') 

n1,o1,d1=from_header('input','n1','o1','d1')
n2,o2,d2=from_header('input','n2','o2','d2')

n1=int(n1)
n2=int(n2)

dat=read('input',(n2,n1))
dat=np.transpose(dat)

# compute data stat
datRMS = np.linalg.norm(dat)/math.sqrt(dat.size)
datMIN = np.min(dat)
datMAX = np.max(dat)

# read parameters
xmin=searchArgv('xmin')
if xmin!=False:
    xmin=float(xmin)
else:
    xmin=o2

xmax=searchArgv('xmax')
if xmax!=False:
    xmax=float(xmax)
else:
    xmax=o2+d2*n2

zmin=searchArgv('zmin')
if zmin!=False:
    zmin=float(zmin)
else:
    zmin=o1

zmax=searchArgv('zmax')
if zmax!=False:
    zmax=float(zmax)
else:
    zmax=o1+d1*n1

ixmin=int((xmin-o2)/d2)
ixmax=int((xmax-o2)/d2)
izmin=int((zmin-o1)/d1)
izmax=int((zmax-o1)/d1)

ixmin=max(0,ixmin)
ixmax=min(n2,ixmax)
izmin=max(0,izmin)
izmax=min(n1,izmax)

dat=dat[izmin:izmax+1, ixmin:ixmax+1]
dat=np.flip(dat,axis=0)

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
    zreverse=1
elif zreverse=='0':
    zreverse=0
else:
    zreverse=1

colormap=searchArgv('colormap')
if colormap!=False:
    colormap=str(colormap)
else:
    colormap='Greys'

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

interpolation=searchArgv('interpolation')
if interpolation==False:
    interpolation='none'

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

colorbar=searchArgv('colorbar')
if colorbar==False:
    colorbar=1
elif colorbar=='0':
    colorbar=0
else:
    colorbar=1

scientific=searchArgv('scientific')
if scientific==False:
    scientific=0
elif scientific=='0':
    scientific=0
else:
    scientific=1

digits=searchArgv('digits')
if digits==False:
    digits=2
else:
    digits=int(digits)

fig, ax = plt.subplots(figsize=(xsize, zsize))

im=ax.imshow(dat,interpolation=interpolation,aspect="auto",extent=[xmin,xmax,zmin,zmax],vmin=min/gain,vmax=max/gain,cmap=colormap)

plt.title(title)
plt.xlabel(xlabel)
plt.ylabel(zlabel)

if colorbar == 1:
    if scientific == 1:
        plt.colorbar(im,format='%.'+str(digits)+'e')
    else:
        plt.colorbar(im,format='%.'+str(digits)+'f')

if zreverse==1:
    plt.gca().invert_yaxis()

out=searchArgv('output')
if out==False: 
    plt.show()
else:
    plt.savefig(out,bbox_inches='tight',format=format)