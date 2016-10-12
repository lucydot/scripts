#!/usr/bin/env python

import argparse
import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import constants
from IPython import embed
from tabulate import tabulate

def read_data(file):
    data =  np.genfromtxt(file, unpack=False, missing_values='nan', 
                          names=True, dtype=None, skip_header=2)
    
    return data

def mandata(data):
  #  for i in range(0,len(data)):
  #      data[i][1]=np.add(data[i][1],0.00626)
  #      data[i][2]=np.add(data[i][2],0.00494)
  #      data[i][3]=np.add(data[i][3], 0.00118)
#    for i in range(0,len(data)):
#        data[i][1]=np.add(data[i][1],0.0177)
#        data[i][2]=np.add(data[i][2],0.0167)
#        data[i][3]=np.add(data[i][3], 0.0)
    return data

def latex_table(data):
    print ('\n')    
    table = [data[i] for i in np.arange(data.shape[0])] 
    headers = data.dtype.names
    latex_table = tabulate(table, headers, tablefmt="latex")
    print (latex_table)
    print ('\n')

def fitdata(xdata, ydata, fitting, symmetry):    
    
    # Symmeterize data so we don't do a ski-jump
    if symmetry:
        negative_xdata=[-value for value in xdata]
        xdata = np.concatenate((xdata,negative_xdata))
        ydata = np.concatenate((ydata,ydata))
        print ("symmeterising data")

    xfitting=[]
    yfitting=[]
    label=[]
    for degree in fitting:
        xfit=np.linspace(xdata[0],xdata[-1],1000, endpoint=True)
        xfitting.append(xfit) 
        yfit=np.polyval(np.polyfit(xdata[:30], ydata[:30], degree), xfit)
        print (xdata[:])
        yfitting.append(yfit)
        print ( "Coeff+Resid for polyfit degree " + str(degree) + " are: " 
                + str(np.polyfit(xdata, ydata, degree, full=True)))
        label.append('fit degree '+ str(degree))
    print ('\n')   
    return xfitting, yfitting, label


def plot(data, xlabel, ylabel, title, xcols, ycols, fitting, 
         nolegend, markers, filename, scatter, symmetry, xlim, ylim):
    
    ## Style
    mpl.rc('font', **{ 'family' : 'serif', 'size' : 8, 
           'serif' : 'Times New Roman' })
    mpl.rc('lines', **{ 'linewidth' : 0.5 })
    
    # to cycle through: mpl.rcParams['axes.color_cycle'] = 
    #['#f79321', '#285fa7', '#215D74']`
    
    # for color maps see
    # http://matplotlib.org/examples/color/colormaps_reference.html
    colour_list = plt.cm.Dark2(np.linspace(0,1,3))
   
    
    plt.figure(figsize = (8.6 / 2.54, 6.0 / 2.54))
    
    if markers is True or scatter is True:
        markers = ["x", "x", "x"]
    else:
        markers = ["",""]
    
    # Plots
    for index,value in enumerate(xcols):
        xname=data.dtype.names[value]
        xdata=data[xname]
        yname=data.dtype.names[ycols[index-1]]
        ydata=data[yname]
        
        # My DIY marker cycler
        if index < len(markers):
            remainder = index
        else:
            remainder = index % len(markers) % index
        
        # My DIY color cycler
        if index < len(colour_list):
            remainder2 = index
        else:
            remainder2 = index % len(colour_list)
        
        if scatter is True:
            plt.scatter(xdata, ydata, marker="x", 
                        label =yname, s=10, color = colour_list[remainder2])
             #marker=markers[remainder],
        else:
            plt.plot(xdata, ydata, marker = markers[remainder], markersize=1,
                     color=colour_list[remainder2], label=yname)
            # marker, markersize, markerfacecolor, markeredgecolor
        
        if fitting is not 'None':
            print("yname: " + yname)
            xfitting, yfitting, label = fitdata(xdata, ydata, fitting, symmetry)
            [plt.plot(xfitting[i],yfitting[i], marker="", linestyle='--', 
            color=colour_list[1]) for i in np.arange(0,len(xfitting))] 
           # ,label=yname+", "+label[i], 
    
    if nolegend is not True:
        plt.legend(loc=2, fontsize=4, ncol=1)
        # bbox_to_anchor=(0.05, 1.05)
    
    for spine in plt.gca().spines.values():
            spine.set_linewidth(0.5)
    # Axes
    if ylim is not None:
        plt.ylim(ylim[0],ylim[1])
    if xlim is not None:
        plt.xlim(xlim[0],xlim[1])
    plt.xticks(np.arange(0.0, max(xdata), 0.01))
    
    #Labels
    if ylabel is not None:
        plt.ylabel(ylabel, fontweight = 'bold', size=12)
    if xlabel is not None:
        plt.xlabel(xlabel, fontweight = 'bold', size=12)
    if title is not None:
        plt.title(title, fontweight = 'bold', fontsize=12)
       
    # File and print
    plt.tight_layout()
    plt.savefig(filename+'.png', format = 'png', dpi = 300)

def main(plots, xlabel, ylabel,
         title, xcols, ycols, fitting, nolegend, markers, filename, scatter,
         symmetry, xlim, ylim):
    
    data=read_data(filename)
    data=mandata(data)    
    # if no columns are specified assume first column is xdata and rest ydata
    if xcols is None:
        xcols=np.zeros(len(data[0])-1, dtype=int)
        ycols=np.arange(1,len(data[0]))

    if 'all' in plots or 'table' in plots:
        latex_table(data)

    if 'all' in plots or 'plot' in plots:
        plot(data, xlabel, ylabel, title, xcols, ycols, fitting, nolegend,
             markers, filename, scatter, symmetry, xlim, ylim)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Generate plots and latex
                                    tables using matplotlib and tabulate""")
    parser.add_argument("-f", "--filename", type=str, default="data",
                        help="filename to read data from")
    parser.add_argument("-p", "--plots", type=str, default="all",
                        help="Plots to generate (as space-separated strings)."
                        "Options are: all, table, plot")
    parser.add_argument("-x", "--xlabel", type=str, default=None,
                        help="Label for x axis.") 
    parser.add_argument("-y", "--ylabel", type=str, default=None,
                        help="Label for y axis.") 
    parser.add_argument("-t", "--title", type=str, default=None,
                        help="Title for graph.")
    parser.add_argument("-xc", "--xcolumns", type=int, default=None, nargs='+',
                        help="Columns of file to use as x values in plot."
                        " Indexing starts at 0." )
    parser.add_argument("-yc", "--ycolumns", type=int, default=None, nargs='+',
                        help="Columns of file to use as y values in plot."
                        " Each must corrspond to a x-column."
                        " Indexing starts at 0." )
    parser.add_argument("-nl", "--nolegend", action="store_true",
                        help=" If selected no legend will be printed.")
    parser.add_argument("-m", "--markers", action="store_true",
                        help=" If selected markers will be printed.")
    parser.add_argument("-fit", "--fitting", type=str, default=None,
                        help="Plot a fit to your data using polyfit"
                        "to polynomial degree specified.")
    parser.add_argument("-s", "--scatter", action='store_true',
                        help=" If selected scatter plot will be made.")
    parser.add_argument("-sym", "--symmetry", action='store_true',
                        help=" If selected data will be symmeterised for any"
                        " fitting. ")
    parser.add_argument("-xlim", "--xlimit", type=float, default=None, nargs='+',
                        help=" x limit for plot")
    parser.add_argument("-ylim", "--ylimit", type=float, default=None, nargs='+',
                        help=" y limit for plot")

    args = parser.parse_args()
    plots = args.plots.split()
    if args.fitting is not None:
        fitting = [int(item) for item in args.fitting.split(',')]
    else:
        fitting = []
    
    main(plots, args.xlabel, args.ylabel, args.title, 
            args.xcolumns, args.ycolumns, fitting, args.nolegend, 
            args.markers, args.filename, args.scatter, args.symmetry, 
            args.xlimit, args.ylimit)
