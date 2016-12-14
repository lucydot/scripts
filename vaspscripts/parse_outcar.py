#!/usr/local/bin/python3
# Simple script for extracting key information from vasp output files

import re

def total_energy(filein):
    # reads from OUTCAR
    read_in = open(filein,'r').read()
    match = re.search("free\s\senergy\s\s*TOTEN\s*=\s*([-.\d\s]+)",read_in )
    if match:
        energy = float(match.groups()[0])
    else:
        energy = 'NaN'
    return energy

def pressure(filein):
    # reads from OUTCAR
    read_in = open(filein,'r').read()
    match = re.search("pressure\s*=\s*([-.\d\s]+)", read_in)
    if match:
        pressure = float(match.groups()[0])
    else:
        pressure = 'NaN'
    return pressure

def total_time(filein):
    # reads from OUTCAR
    read_in = open(filein,'r').read()
    match = re.search("time\sused\s\(sec\):\s*([-.\d\s]+)",read_in)
    if match:
        time = float(match.groups()[0])
    else:
        time = 'NaN'
    return time

def ave_SCF_time(filein, skipSCF=0):
    # reads from OUTCAR
    # Adapted from github:JMSKelton
    # skipSCF may be useful when doing hybrid of algo=fast calculations
    # where the first five cycles are of a different length
    
    scf_times=[]

    with open(filein, 'r') as read_in:
        for line in read_in:
            match = re.compile("LOOP\:\s+cpu time\s+\d+\.\d+\:\s+real time\s+(?P<t_scf>\d+\.\d+)").search(line)

            if match:
                scf_times.append(float(match.group('t_scf')))

    if skipSCF > 0:
        if len(scf_times) > skipSCF:
            scf_times = scf_times[skipSCF:]
        else:
            print ("WARNING: _average_time: Number of SCF steps {0} <= skip SCF {1}".format(len(scf_times), skipSCF))

    if len(scf_times) > 0:
        ave_SCF_time = round(sum(scf_times) / len(scf_times))
    
    else:
        ave_SCF_time = 'Error'
    return (ave_SCF_time)

def SCF_steps(filein):
    scf_times=[]
    with open(filein, 'r') as read_in:
        for line in read_in:
            match = re.compile("LOOP\:\s+cpu time\s+\d+\.\d+\:\s+real time\s+(?P<t_scf>\d+\.\d+)").search(line)

            if match:
                scf_times.append(float(match.group('t_scf')))

    return len(scf_times)

def kpar(filein):
    # reads from OUTCAR
    return 0

def ncore(filein):
    # reads from OUTCAR
    return 0

def cores(filein):
    # reads from OUTCAR
    return 0

def nbands(filein):
    # reads from OUTCAR
    return 0

def nkpts(filein):
    # reads from OUTCAR
    # KPOINTS calculated after symmetry considerations. May not be product
    # of kmesh dimensions
    return 0

def nion(filein, atoms_in_unit_cell = None):
    # reads from OUTCAR
    if atoms_in_unit_cell:
        # calculate how many unit cells there are
    return 0

def nelect(filein):
    # reads from OUTCAR
    return 0

def kmesh(filein):
    # reads from ....... KPOINTS!


