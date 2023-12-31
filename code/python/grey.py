#!/usr/bin/env python

import matplotlib
#matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from seppyio import *
from matplotlib.patches import FancyBboxPatch 
from mpl_toolkits.axes_grid1 import make_axes_locatable

filename=searchArgv('hfile') 

if searchArgv('o1')==False:
    n1,o1,d1=from_header('hfile','n1','o1','d1')
    n2,o2,d2=from_header('hfile','n2','o2','d2')
else:
    o1,d1=get_param('o1','d1')
    n1,n2=from_header('hfile','n1','n2')
    o2,d2=get_param('o2','d2')

n1=int(n1)
n2=int(n2)

a=read('hfile',(n2,n1))

temp=searchArgv('scalefactor')
if temp!=False:
    temp=float(temp)
    a=a*temp

min1=searchArgv('min1')
if min1!=False:
    min1=float(min1)
    b1=int((min1-o1)/d1)
else:
    b1=0

max1=searchArgv('max1')
if max1!=False:
    max1=float(max1)
    e1=int((max1-o1)/d1)
    o1=min1
else:
    e1=n1
n1=e1-b1

min2=searchArgv('min2')
if min2!=False:
    min2=float(min2)
    b2=int((min2-o2)/d2)
else:
    b2=0

max2=searchArgv('max2')
if max2!=False:
    max2=float(max2)
    e2=int((max2-o2)/d2)
    o2=min2
else:
    e2=n2
n2=e2-b2

a=a[b2:e2,b1:e1]

transp=searchArgv('transp')
if transp=='y': 
    a=np.transpose(a)
    n=n1
    n1=n2
    n2=n
    o=o1
    o1=o2
    o2=o
    d=d1
    d1=d2
    d2=d

temp=searchArgv('maxval')
if temp==False: 
    maxval=np.max(a) 
else:
    maxval=float(temp)

temp=searchArgv('minval')
if temp==False: 
    minval=np.min(a) 
else:
    minval=float(temp)

temp=searchArgv('colormap')
if temp==False: 
    colormap='gray' 
else:
    colormap=temp

temp=searchArgv('aspect')
if temp==False:
    ar='auto'
else:
    ar=float(temp)

w,h=get_param("figwidth","figheight")
if w!=False:
    fig, ax = plt.subplots(ncols=1,figsize=(w,h))
else:
    fig, ax = plt.subplots(ncols=1)

im=ax.imshow(a,interpolation='none',aspect="auto",extent=[o1,o1+(n1)*d1,o2+(n2)*d2,o2],vmin=minval,vmax=maxval,cmap=colormap)

temp=searchArgv('fontsize')
if temp==False:
    sz=15
else:
    sz=float(temp)

ax.tick_params(labelsize=sz,pad=0.5*sz)
#ax.get_xaxis().set_ticks([])
#ax.get_yaxis().set_ticks([])
#plt.setp(ax.get_xticklabels(), visible=False)
#plt.setp(ax.get_yticklabels(), visible=False)

temp=searchArgv('xlabel')
if temp!=False: 
    xlabel=temp
    ax.set_xlabel(xlabel,fontsize=sz)

temp=searchArgv('ylabel')
if temp!=False: 
    ylabel=temp
    ax.set_ylabel(ylabel,fontsize=sz)

temp=searchArgv('title')
if temp==False: 
    title=filename
else:
    title=temp
ax.set_title(title,fontsize=sz)

temp=searchArgv('colorbar')
if temp=='y':
    divider = make_axes_locatable(ax)
    barloc=searchArgv("barloc")
    if barloc==False:
        barloc="right"
    barsize=searchArgv("barsize")
    if barsize==False:
        barsize="3%"
    barpad=get_param("barpad")
    if barpad==False:
        barpad=0.2
    cax = divider.append_axes(barloc,size=barsize,pad=barpad)
    ave=0.5*(minval+maxval)
    cbar=fig.colorbar(im,cax=cax,ticks=[minval,0.5*(minval+ave),ave,0.5*(ave+maxval),maxval])
    cbar.ax.tick_params(labelsize=sz,pad=0.5*sz) 
    temp=searchArgv('barlabel')
    if temp!=False:
        cbar.ax.set_ylabel(temp,fontsize=sz)

annotate=searchArgv('annotate')
if annotate=='box':
    lowerleft=get_array("lowerleft")
    boxSize=get_array("boxSize")
    boxPad=searchArgv("boxPad")
    boxColor=searchArgv("boxColor")
    boxLineWidth=searchArgv("boxLineWidth")
    p_bbox=FancyBboxPatch(lowerleft,boxSize[0],boxSize[1],ec="k",fc="none",linewidth=boxLineWidth)
    p_bbox.set_boxstyle("round",pad=boxPad);
    ax.add_patch(p_bbox)

if annotate=='arrow':
    tipx=get_array('arrowtipx')
    tipy=get_array('arrowtipy')
    textx=get_array('arrowtextx')
    texty=get_array('arrowtexty')
    narrow=len(tipx)
    arrowColor=get_sarray('arrowColor')
    if arrowColor==False:
        arrowColor=[]
        for i in range(narrow):
            arrowColor.append('k')
    arrowTitle=get_sarray('arrowTitle')
    if arrowTitle==False:
        arrowTitle=[]
        for i in range(narrow):
            arrowTitle.append('')
    for i in range(narrow):
        ax.annotate(arrowTitle[i],xy=(tipx[i],tipy[i]),xytext=(textx[i],texty[i]),arrowprops=dict(arrowstyle='->',facecolor=arrowColor[i]))

grid=searchArgv('grid')
if grid=='y':
    ax.grid()

gridx=searchArgv('gridx')
if gridx=='y':
    ax.grid(axis='x')

gridy=searchArgv('gridy')
if gridy=='y':
    ax.grid(axis='y')

temp=searchArgv('out')
if temp==False: 
    plt.show()
else:
    plt.savefig(temp,bbox_inches='tight')

