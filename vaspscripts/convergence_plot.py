#!/usr/local/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import re
import argparse
from IPython import embed

def total_energy():
    free_energies = main(begin=1,save='False', kind='electronic',filename='OUTCAR')
    return free_energies[-1]

def set_axlims(array, marginfactor):
    """
        Fix for a scaling issue with matplotlibs scatterplot and small values.
        Takes in a numpy array, and a marginfactor (float).
        A marginfactor of 0.2 would for example set a 20% border distance on both sides.
        Output:[bottom,top]
        To be used with .set_ylim(bottom,top)
    """
    minv = array.min()
    maxv = array.max()
    datarange = maxv-minv
    border = abs(datarange*marginfactor)
    maxlim = maxv+border
    minlim = minv-border

    return minlim,maxlim

def main(begin=1, save=False, kind='electronic',filename='OUTCAR'):

    read_in = open(filename,'r').read()

    if kind == 'electronic':
        free_energies = re.findall("free\senergy\s*TOTEN\s*=\s*([-.\d\s]+)",read_in)     
    if kind == 'ionic':
        free_energies = re.findall("free\s\senergy\s\s*TOTEN\s*=\s*([-.\d\s]+)",read_in) 
    
    free_energies = [i for i in free_energies if i != ' ']
    free_energies = [float(i) for i in free_energies][begin-1:]
    print (free_energies)

    minlim, maxlim = set_axlims(np.array(free_energies), 0.05)

    x_coords = range(1, len(free_energies)+1)
    plt.scatter(x_coords,free_energies)
    
    axes = plt.gca()
    axes.set_ylim([minlim,maxlim])
    plt.savefig('energies.png') 
   
    plt.clf()

    energies_difference = np.diff(free_energies)
    print (energies_difference)

    minlim, maxlim = set_axlims(np.array(energies_difference), 0.05)
    x_coords = range(1,len(energies_difference)+1)
    plt.scatter(x_coords,energies_difference)
    
    axes = plt.gca()
    axes.set_ylim([minlim,maxlim])

    plt.savefig('convergence.png')
    
    return free_energies

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot energy convergence')
    parser.add_argument('-t', '--type', default='electronic', help='Select electronic or ionic convergence (default:electronic)')
    parser.add_argument('-b', '--begin', default=1, help='step to start plot at (default:1)', type=int)
    parser.add_argument('-s', '--save', default=False, help='save plot (default:False)', type=bool)
    parser.add_argument('-f','--filename', default='OUTCAR', help='File to read from')
    args = vars(parser.parse_args())
    main(begin=args["begin"],save=args["save"],kind=args["type"], filename=args["filename"])
   

