#!/usr/local/bin/python3
import re
import sh
from IPython import embed

data_folder="/Users/lucydot/data/151/"
calculation_folders=("111-444","111-555","111-666","111-777","111-888","111-999","111-101010","111-111111","111-121212","111-131313","111-141414","111-151515")
results_folder="/Users/lucydot/projects/nr-recombo/results/convergence_tests/"
filename="MAPI_kpt_scaling"

with open(results_folder+filename,"a") as myfile:
    myfile.write("calculation energy pressure time \n")
    for calculation in calculation_folders:
        sh.cd(data_folder+calculation)
        read_in = open('OUTCAR','r').read()
        energy = re.search("free\s\senergy\s\s*TOTEN\s*=\s*([-.\d\s]+)",read_in )
        pressure = re.search("pressure\s*=\s*([-.\d\s]+)", read_in)
        time = re.search("time\sused\s\(sec\):\s*([-.\d\s]+)",read_in)
        myfile.write(calculation + " " + (energy.groups()[0]) + " " + (pressure.groups()[0]) + " " + (time.groups()[0]) + "\n" ) 



