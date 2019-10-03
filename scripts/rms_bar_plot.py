#!/usr/bin/env python

# This script is used for plotting the orbit RMS from precise orbit determiantion.
# Please change the path of file and file name. The figure is saved as *.png.

# Output file name format is expressed as year+DOY+hour(system time).

# Authors: Tzupang Tseng, Carl Wang and Salim Masoumi
#
# Date: 11-06-2019
#
#

import os
import sys
import argparse
from argparse import ArgumentParser
import time
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import gc
import math

# File-related statements
os.chdir('/data/test/')

# Command line argument
parser = argparse.ArgumentParser()

parser.add_argument('-i', type=str, required=True, help='Input file name')

args = parser.parse_args()

inputfile = args.i

# End of command line argument


infile = open(inputfile, 'r')
mat_1 = np.loadtxt(infile)
year = time.strftime("%Y")
doy  = time.strftime("%j")
hr   = time.strftime("%H")
minu = time.strftime("%M")
# Output file name format is expressed as year+DOY+hour(system time).
file_plt_1='Orbrms'+year+doy+'_'+hr+'.png'
########################################################

#-------plot 1-1------------------------------------------------------------
prn = mat_1[:,1]
plt.figure(figsize=(13,8))
plt.subplot(111)
plt.title("Satellite orbit residuals",fontsize=20)
bar_width = 0.25
for I in range(1,int(max(prn))+1):
    r = []
    t = []
    n = []
    for x in range (0,len(prn)):
        if prn[x] == I:
           r.append(mat_1[x,10]*100)
           t.append(mat_1[x,11]*100)
           n.append(mat_1[x,12]*100)
    r_mean= np.mean(r,dtype=np.float64)
    t_mean= np.mean(t,dtype=np.float64)
    n_mean= np.mean(n,dtype=np.float64)
    r_std = np.std(r,dtype=np.float64)
    t_std = np.std(t,dtype=np.float64)
    n_std = np.std(n,dtype=np.float64)
    r_rms = math.sqrt(r_mean**2+r_std**2)
    t_rms = math.sqrt(t_mean**2+t_std**2)
    n_rms = math.sqrt(n_mean**2+n_std**2)

    if I ==1:
       plt.bar(I-bar_width, r_rms, bar_width,alpha=0.8, color='r', align='center', label='Radial')
       plt.bar(I, t_rms, bar_width, alpha=0.8, color='g', align='center', label='Along-track')
       plt.bar(I+bar_width, n_rms, bar_width,alpha=0.8, color='b', align='center', label= 'Cross-track')
       plt.legend()
    else:
        plt.bar(I-bar_width, r_rms, bar_width, color='r', align='center')
        plt.bar(I, t_rms, bar_width, color='g', align='center')
        plt.bar(I+bar_width, n_rms, bar_width, color='b', align='center')


plt.xlim(0, int(max(prn))+1)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.ylabel('RMS (cm)',fontsize=20)
plt.xlabel('Satellite PRN number',fontsize=20)
plt.title('Satellite orbit RMS',fontsize=20)
plt.grid()
plt.savefig(file_plt_1,dpi=150)

#plt.show()

#plt.clf()

