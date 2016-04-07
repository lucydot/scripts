#!/usr/local/bin/python
import matplotlib.pyplot as plt
import re
import argparse

parser = argparse.ArgumentParser(description='Plot energy convergence')
parser.add_argument('-t', '--type', default='ionic', help='Select electronic or ionic convergence (default:ionic)')
args = vars(parser.parse_args())

read_in = open('OUTCAR','r').read()

if args['type'] == 'electronic':
    free_energies = re.findall("free\senergy\s*TOTEN\s*=\s*([-.\d\s]+)",read_in) # 2 spaces for ionic energy step, 1 space if wanted electronic steps
if args['type'] == 'ionic':
    free_energies = re.findall("free\s\senergy\s\s*TOTEN\s*=\s*([-.\d\s]+)",read_in) # 2 spaces for ionic energy step, 1 space if wanted electronic steps


[float(i) for i in free_energies]

x_coords = range(1, len(free_energies)+1)

print free_energies

sc = plt.scatter(x_coords,free_energies)
plt.show() 
