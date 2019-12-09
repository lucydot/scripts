#! /usr/bin/env python

import numpy
import matplotlib.pyplot as plt
from pylab import *
import math
from optparse import OptionParser


def file_len(fname):
    """ Returns the number of lines in a file """
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

parser = OptionParser()
parser.add_option("-x", "--xaxis",
                  action="store", type="string", dest="x_lab", default="X [arb.]",
                  help="x axis label [default: X [arb.]]")

parser.add_option("-y", "--yaxis",
                  action="store", type="string", dest="y_lab", default="Y [arb.]",
                  help="y axis label [default: Y [arb.]]")

parser.add_option("-f", "--file",
                  action="store", type="string", dest="file", default="Bands.dat",
                  help="Path to input file [default: Phonons.dat]")

parser.add_option("-n", "--number",
                  action="store", type="int", dest="n_data", default="1",
                  help="Number of datasets to plot [default = 1]")

parser.add_option("-i", "--fill",
                  action="store", type="string", dest="fill_in", default="false",
                  help="Fill in unber the lines (true/false) [default = false]")

parser.add_option("-t", "--title",
                  action="store", type="string", dest="title", default="Band diagram",
                  help="Title [default: Band diagram]")
parser.add_option("-s", "--save",
                  action="store", type="string", dest="filename",
                  help="Save as filename")
# Add further options here
(options, args) = parser.parse_args()
# Ensure that the input files exist

try:
 with open(options.file) as tmp: pass
except IOError as e:
 print 'Could not find a coordinates file, this file is specified with the -f flag, default Data.dat'

# Read the data
Dataset = np.zeros(shape=(file_len(options.file),options.n_data+1))
f = open(options.file,"r")
lines = f.readlines()
f.close()
i = 0
for line in lines:
    inp = line.split()
    Dataset[i,:] = inp[:]
    i = i + 1
#Plotting
#STYLE
xticklines = getp(gca(), 'xticklines')
yticklines = getp(gca(), 'yticklines')
xgridlines = getp(gca(), 'xgridlines')
ygridlines = getp(gca(), 'ygridlines')
xticklabels = getp(gca(), 'xticklabels')
ygridlines = getp(gca(), 'ygridlines')
xticklabels = getp(gca(), 'xticklabels')
yticklabels = getp(gca(), 'yticklabels')
setp(xticklines, 'linewidth', 3)
setp(yticklines, 'linewidth', 3)
#setp(xgridlines, 'linestyle', '-')
#setp(ygridlines, 'linestyle', '-')
setp(yticklabels, 'color', 'Black', fontsize='medium')
setp(xticklabels, 'color', 'Black', fontsize='medium')
plt.rc('axes', color_cycle=['FF8000', 'CCCC00', '3399FF', '99004C'])
# AXIS TITLES
plt.xlabel(options.x_lab,fontsize=18)
plt.ylabel(options.y_lab,fontsize=18)
# TITLE
plt.title(options.title)
# Plot 
i = 1
while i <= options.n_data :
    plt.plot(Dataset[:,0],Dataset[:,i],linewidth=1.5)
    if options.fill_in == "true":
 	plt.fill(Dataset[:,0],Dataset[:,i],'#FF8000',alpha=0.5)

    i = i + 1
plt.grid(True)

if options.filename is not None:
    plt.savefig(options.filename)

else:
    plt.show()


