#!/usr/local/bin/python
import matplotlib.pyplot as plt
import re
import argparse

parser = argparse.ArgumentParser(description='Plot energy convergence')
parser.add_argument('-t', '--type', default='electronic', help='Select electronic or ionic convergence (default:electronic)')
parser.add_argument('-b', '--begin', default=1, help='step to start plot at (default:1)', type=int)
parser.add_argument('-s', '--save', default=False, help='save plot (default:False)', type=bool)
args = vars(parser.parse_args())

read_in = open('OUTCAR','r').read()

if args['type'] == 'electronic':
    free_energies = re.findall("free\senergy\s*TOTEN\s*=\s*([-.\d\s]+)",read_in) # 2 spaces for ionic energy step, 1 space if wanted electronic steps
if args['type'] == 'ionic':
    free_energies = re.findall("free\s\senergy\s\s*TOTEN\s*=\s*([-.\d\s]+)",read_in) # 2 spaces for ionic energy step, 1 space if wanted electronic steps


[float(i) for i in free_energies]

x_coords = range(1, len(free_energies)+1)

print free_energies

plt.scatter(x_coords[args['begin']-1:],free_energies[args['begin']-1:])
if args['save'] == True:
    plt.savefig('convergence.png') 
plt.show()

