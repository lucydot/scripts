#!/usr/local/bin/python
import matplotlib.pyplot as plt
import re
import argparse

def total_energy():
    plt, free_energies = main(begin=1,save='False', kind='electronic',file='OUTCAR')
    return free_energies[-1]

def main(begin=1, save=False, kind='electronic',file='OUTCAR'):

    read_in = open(file,'r').read()

    if kind == 'electronic':
        free_energies = re.findall("free\senergy\s*TOTEN\s*=\s*([-.\d\s]+)",read_in)     
    if kind == 'ionic':
        free_energies = re.findall("free\s\senergy\s\s*TOTEN\s*=\s*([-.\d\s]+)",read_in) 

    [float(i) for i in free_energies]

    x_coords = range(1, len(free_energies)+1)

    plt.scatter(x_coords[begin-1:],free_energies[begin-1:])
    
    if kind == True:
        plt.savefig('convergence.png') 
    
    return plt, free_energies

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot energy convergence')
    parser.add_argument('-t', '--type', default='electronic', help='Select electronic or ionic convergence (default:electronic)')
    parser.add_argument('-b', '--begin', default=1, help='step to start plot at (default:1)', type=int)
    parser.add_argument('-s', '--save', default=False, help='save plot (default:False)', type=bool)
    args = vars(parser.parse_args())
    plt, free_energies = main(begin=args["begin"],save=args["save"],kind=args["type"], file="OUTCAR")
    print (free_energies)
    plt.show()
   

