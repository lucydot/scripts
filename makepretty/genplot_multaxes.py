#!/usr/bin/env python
# Beginnings of multiple reading in.

import argparse
import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import constants
from IPython import embed
from tabulate import tabulate

def read_data(files):
    data=[]
    for index,file in enumerate(files):
        data.append(np.genfromtxt(file, unpack=False, missing_values='nan', 
                          names=True, dtype=None, skip_header=2))
    
    return data

def mandata(data):
    # do stuff to data
    return data


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


def plot(data, args):

    ## Style
    fontsize=3
    mpl.rc('font', **{ 'family' : 'serif', 'size' : fontsize, 
           'serif' : 'Times New Roman' })
    mpl.rc('lines', **{ 'linewidth' : 0.5 })
    mpl.rcParams['axes.linewidth'] = 0.5
    mpl.rcParams['xtick.major.width'] = 0.5
    mpl.rcParams['ytick.major.width'] = 0.5
    # changing things globally...can also set for individual axes.
    
    # to cycle through: mpl.rcParams['axes.color_cycle'] = 
    #['#f79321', '#285fa7', '#215D74']`
    
    # for color maps see
    # http://matplotlib.org/examples/color/colormaps_reference.html
    color_list = plt.cm.Dark2(np.linspace(0,1,3))
    
    if args["markers"] is True or args["scatter"] is True:
        marker_list = ["x", "x", "x"]
    else:
        marker_list = ["","","",""]
    
    # need to input tiling
    f, axarray = plt.subplots(2,2,figsize = (6 / 2.54, 6 / 2.54),
                              sharey = True, sharex=True)
    # see http://matplotlib.org/examples/pylab_examples/subplots_demo.html
    # sharey = row for example
    
    # need to input tiling positions
    order = [(0,0),(0,1),(1,0),(1,1)]
    
    labels = ['PBEsol','PBEsol + SoC','HSE06']
    
   # Plot it up
    for series,dat in enumerate(data):
        for index,value in enumerate(args["xcolumns"]):
            xname=data[series].dtype.names[value]
            xdata=data[series][xname]
            yname=data[series].dtype.names[args["ycolumns"][index]]
            ydata=data[series][yname]
            if args["manual_labels"] is True:
                name = labels[index]
            else:
                name=yname
      
            # Fitting. Index labels dataset, i labels fitting.
            # If you want to add labelling for diff fits set label=label[i]
            if args["fitting"] is not []:
                print("yname: " + yname)
                xfitting, yfitting, label = fitdata(xdata, ydata, args["fitting"], args["symmetry"])
                [axarray[order[series]].plot(xfitting[i],yfitting[i], marker="", label=None,
                color=color_list[index]) 
                for i in np.arange(0, len(xfitting))  ] 
            
            if args["scatter"] is True:
                axarray[order[series]].scatter(xdata, ydata, marker=marker_list[index],
                            color = color_list[index], label=label)
            else:
                axarray[order[series]].plot(xdata, ydata, marker = marker_list[index], markersize=1,
                         color=color_list[index], label=name)
                # marker, markersize, markerfacecolor, markeredgecolor
     
        axarray[order[series]].legend(frameon=False,loc=2, fontsize=fontsize, ncol=1)
        
       # bbox_to_anchor=(0.05, 1.05)
            
    # Axes
    if args["ylimit"] is not None:
        axarray[series].ylim(args["ylimit"][0],args["ylimit"][1])
    if args["xlimit"] is not None:
        axarray[series].xlim(args["xlimit"][0],args["xlimit"][1])
#    axarray[series].xticks(np.arange(0.0, max(xdata), 0.01))
    
    #Labels
    if args["ylabel"] is not None:
        axarray[series].ylabel(args["ylabel"][series], size=fontsize)
    if args["xlabel"] is not None:
        axarray[series].xlabel(args["xlabel"][series], size=fontsize)
    if args["title"] is not None:
        axarray[series].title(args["title"][series], fontsize=fontsize)
    
    if args["plotname"] is not None:
        plotname = args["plotname"]
    else:
        plotname=args["filename"][0]

    # File and print
    plt.tight_layout()
    plt.savefig(plotname+'.png', format = 'png', dpi = 300)
    plt.close("all")

def main(args):
    
    data=read_data(args["filename"])
    data=mandata(data)    
    
    # if no columns are specified assume first column is xdata and rest ydata
    if args['xcolumns'] is None:
        args['xcolumns'] = np.zeros(len(data[0][0])-1, dtype=int)
        args['ycolumns'] = np.arange(1,len(data[0][0]))

    plot(data, args)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Generate plots and latex
                                    tables using matplotlib and tabulate""")
    parser.add_argument("-f", "--filename", type=str, default="data",nargs='+',
                        help="filenameis to read data from")
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
    parser.add_argument("-fit", "--fitting", type=str, default=[],
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
    parser.add_argument("-ml", "--manual_labels", action='store_true',
                        help=" If selected, labels can be input manually")
    parser.add_argument("-pn", "--plotname", type=str, default=None,
                        help="plotname for saving")
    args = vars(parser.parse_args())
    
    main(args) 
