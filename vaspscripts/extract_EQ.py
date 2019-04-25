#!/usr/bin/env python3

from pymatgen.io.vasp.inputs import Poscar
from scipy.linalg import solve
from io import StringIO
from os import remove
from argparse import ArgumentParser
from re import findall
from csv import writer
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


##### functions read_poscar and get_init_fin written by Sunghyun Kim ######

def read_poscar(i_path, l_get_sorted_symbols=False):
    poscar = Poscar.from_file("{}".format(i_path))
    struct = poscar.structure
    if l_get_sorted_symbols:
        return struct, poscar.site_symbols
    else:
        return struct


def get_init_fin(i_file, f_file):

    # A. Alkauskas, Q. Yan, and C. G. Van de Walle, Physical Review B 90, 27 (2014)

    struct_i, sorted_symbols = read_poscar(i_file, True)
    struct_f, sorted_symbols = read_poscar(f_file, True)
    delta_R = struct_f.frac_coords - struct_i.frac_coords
    delta_R = (delta_R + 0.5) % 1 - 0.5
    
    lattice = struct_i.lattice.matrix
    delta_R = np.dot(delta_R, lattice)


    masses = np.array([spc.atomic_mass for spc in struct_i.species])
    delta_Q2 = masses[:,None] * delta_R ** 2
    return (np.sqrt(delta_Q2.sum()))


def make_all_POSCAR(outcar="OUTCAR",headings_file="POSCAR",fileout="POSCAR_"):
    num_iterations = len(findall("free\s\senergy\s\s*TOTEN\s*=\s*([-.\d\s]+)",open(outcar,'r').read()))
    for i in range (1,num_iterations+1):
        make_POSCAR(i,outcar=outcar,headings_file=headings_file,fileout=fileout)
        
def make_POSCAR(iteration,outcar="OUTCAR",headings_file="POSCAR",fileout="POSCAR_"):
    
    # read in data from outcar and extract list of atomic positions at a given ionic step
    data = open(outcar,'r').read()
    position_data = data.split(str(iteration)+"(",1)[-1].split("total drift,1")[-1].split("total drift")[0].split("Angst)")[-1]
    df = pd.read_csv(StringIO(position_data),delimiter="\n",names="r",skiprows=2)
    positions = [df.r[i].split()[:3] for i in range(df.r.count()-1)]
    position_floats = np.array([[float(value) for value in vector] for vector in positions])
    
    # read the cell data from the POSCAR
    i=0
    cell_matrix = np.zeros(shape=(3,3))
    with open(headings_file,'r') as filein:
        for line in filein:
            if isfloat(line.strip()):
                scaling = float(line.strip())

            elif ([isfloat(item) for item in line.strip().split()] == [True]*len(line.strip().split())) & (i<3):
                cell_vector=[float(value) for value in (line.strip().split())]
                cell_matrix[i]=cell_vector
                i+=1      
    cell_matrix=np.swapaxes(cell_matrix,0,1)
 
    # use cell data to convert positions in angstrom to fractional coordinates
    position_fractional = []
    for position in position_floats:
        position_fractional.append(list(solve(cell_matrix,position)))
    
    # read the POSCAR heading information and copy to new fileout
    with open(headings_file,'r') as filein, open(fileout+str(iteration),'w') as fileout:
        copy = True
        for line in filein:                
            if line.strip() == "Direct":
                fileout.write(line)
                copy == False
                break
            elif copy:
                fileout.write(line)
        # write the atomic positions to fileout
        for position in position_fractional:
            fileout.write(str(position[0])+"     "+str(position[1])+"     "+str(position[2])+"\n")

def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

def find_EQ(outcar="OUTCAR",headings_file="POSCAR",fileout="POSCAR_temp_",keep_poscars=False):
    
    free_energies = [float(x) for x in findall("free\s\senergy\s\s*TOTEN\s*=\s*([-.\d\s]+)",open(outcar,'r').read())]
    num_iterations = len(free_energies)
    
    make_all_POSCAR(outcar=outcar, headings_file=headings_file, fileout="POSCAR_temp_")
    
    Qs=[]
    for i in range(1,num_iterations+1):
        Qs.append(get_init_fin("./POSCAR","./POSCAR_temp_"+str(i)))
    
    if not keep_poscars:
        for i in range(1,num_iterations+1):
            remove("./POSCAR_temp_"+str(i))

    plt.scatter(Qs,free_energies)
    plt.xlabel("Q")
    plt.ylabel("Energy (eV)")
    plt.savefig("E-Q.png")
    
    rows = zip(Qs,free_energies)
    with open("E-Q.txt", "w") as fileout:
        for row in rows:
            writer(fileout).writerow(row)
            
    return free_energies,Qs

if __name__== '__main__':
    parser = ArgumentParser(description='extract atomic positions from OUTCAR at particular iteration and use this to plot energy vs Q (configuration coordinate)')
    parser.add_argument('-fi','--filein', default='OUTCAR',type=str, help='name of file (usually OUTCAR) you would like to extract positions from')
    parser.add_argument('-hi','--heading', default='POSCAR',type=str, help='name of file (usually POSCAR) you would like to extract headings from')
    parser.add_argument('-kp','--keep_poscars',default=False,type=bool,help='if True, the extracted POSCARs will be kept')
    args = parser.parse_args()

    find_EQ(outcar=args.filein,headings_file=args.heading,keep_poscars=args.keep_poscars)

    

