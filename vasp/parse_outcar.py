#simple script for extracting key information from vasp output files

import re

def total_energy(filein):
    # reads from OUTCAR
    read_in = open(filein,'r').read()
    match = re.search("free\s\senergy\s\s*TOTEN\s*=\s*([-.\d\s]+)",read_in )
    if match:
        energy = float(match.groups()[0])
    else:
        energy = None
    return energy

def pressure(filein):
    # reads from OUTCAR
    read_in = open(filein,'r').read()
    match = re.search("pressure\s*=\s*([-.\d\s]+)", read_in)
    if match:
        pressure = float(match.groups()[0])
    else:
        pressure = None
    return pressure

def total_time(filein):
    # reads from OUTCAR
    read_in = open(filein,'r').read()
    match = re.search("time\sused\s\(sec\):\s*([-.\d\s]+)",read_in)
    if match:
        time = float(match.groups()[0])
    else:
        time = None
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
        ave_SCF_time = None
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
    read_in = open(filein,'r').read()
    match = re.search("distrk:\s*each\sk-point\son\s*\d*\scores,\s*([\d]*)",read_in)
    if match:
        kpar = int(match.groups()[0])
    else:
        kpar = None
    return kpar

def ncore(filein):
    read_in = open(filein,'r').read()
    match = re.search("distr:\s*one\sband\s*on\sNCORES_PER_BAND=\s*([\d]*)", read_in)
    if match:
        ncore = int(match.groups()[0])
    else:
        ncore = None
    return ncore

def npar(filein):
    read_in = open(filein,'r').read()
    match = re.search("distr:\s*one\sband\s*on\sNCORES_PER_BAND=\s*\d*\scores,\s*([\d]*)", read_in)
    if match:
        npar = int(match.groups()[0])

    return npar

def cores(filein):
    # reads from OUTCAR
    kpa = kpar(filein)
    read_in = open(filein,'r').read()
    match = re.search("distrk:\s*each\sk-point\son\s*([-.\d]+)",read_in)
    if match:
        kcores = int(match.groups()[0])
    else:
        kcores = None
    if kpa and kcores:
        cores = kpa*kcores
    else:
        cores= None
    return cores

def nbands(filein):
    # reads from OUTCAR
    read_in = open(filein,'r').read()
    match = re.search("NBANDS=\s*([\d]*)",read_in)
    if match:
        nbands = int(match.groups()[0])
    else:
        nbands = None
    return nbands

def nkpts(filein):
    # reads from OUTCAR
    # KPOINTS calculated after symmetry considerations. May not be product
    # of kmesh dimensions
    read_in = open(filein,'r').read()
    match = re.search("NKPTS\s=\s*([\d]*)",read_in)
    if match:
        nkpts = int(match.groups()[0])
    else:
        nktps = None
    return nkpts

def nion(filein):
    # reads from OUTCAR
    read_in = open(filein,'r').read()
    match = re.search("NIONS\s=\s*([\d]*)" ,read_in)
    if match:
        nions = int(match.groups()[0])
    else:
        nions = None
    return nions

def ncells(filein, atoms_in_unit_cell):
    # reads from OUTCAR
    read_in = open(filein,'r').read()
    match = re.search("NIONS\s=\s*([\d]*)" ,read_in)
    if match:
        nions = int(match.groups()[0])
    else:
        nions = None
    if nions:
        unit_cells = float(nions / atoms_in_unit_cell)
    else:
        unit_cells = None
    return unit_cells

def nelect(filein):
    # reads from OUTCAR
    read_in = open(filein,'r').read()
    match = re.search("NELECT\s=\s*([.\d]*)", read_in)
    if match:
        nelect = float(match.groups()[0])
    else:
        nelect = None
    return nelect

def kmesh(filein):
    # reads from ....... KPOINTS!
    read_in = open(filein,'r').read()
    match = re.search("centered\sgrid\s*([.\d\s]*)" ,read_in)
    if match:
        kmesh = [int(i) for i in match.groups()[0].split()]
    else:
        kmesh = None
    
    return kmesh

