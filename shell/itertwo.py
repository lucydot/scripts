#!/usr/local/bin/python3
import re
import sh
from vaspscripts import parse_outcar

root_folder="/Users/lucydot/data/151/"
folders=("111-444","111-555","111-666","111-777","111-888","111-999","111-101010","111-111111","111-121212","111-131313","111-141414","111-151515")
results_file="/Users/lucydot/projects/nr-recombo/results/convergence_tests/XXXX"

with open(results_file,"a") as myfile:
    myfile.write("calculation energy pressure ave_SCF_time time \n")
    for folder in folders:
        sh.cd(root_folder+folder)
        energy = parse_outcar.total_energy('OUTCAR')
        pressure = parse_outcar.pressure('OUTCAR')
        time = parse_outcar.total_time('OUTCAR')
        ave_SCF_time = parse_outcar.ave_SCF_time('OUTCAR')
        myfile.write(folder + " " + str(energy) + " " + str(pressure) + " " + str(ave_SCF_time) + " " +time + "\n") 


