#!/usr/bin/env python

import numpy as np
import os
from IPython import embed
import re

folder = "/Users/lucydot/data/all/149/"
mode = "M"
numbers = (1,2,3,4,5,6,7,10,11,12,13)

VBM=[]
CBM=[]
Pb=[]

for n in numbers:
    zeron = format(n, '03')
    print ("In folder "+zeron)
    os.chdir(folder+zeron)
    os.system("bandenergies")    
    x=np.genfromtxt('OUTCAR - Band Energies.csv',skip_header=2, 
            unpack=False, delimiter='","')
    
    outcar=open('OUTCAR', "r").read()
    Pb_energy=re.search(r"3d\s*([-.\d]*)",outcar)
    Pb_energy=float(Pb_energy.group(1))

    VBM.append(x[199,1])
    CBM.append(x[200,1])
    Pb.append(Pb_energy)


os.chdir(folder)
newfolder="Analysis"
os.system("mkdir "+newfolder)
os.chdir(folder+newfolder)

adjusted_CBM = np.subtract(CBM,CBM[0])
adjusted_Pb = np.subtract(Pb,Pb[0])

Pb_VBM=np.subtract(VBM,adjusted_Pb)
Pb_CBM=np.subtract(adjusted_CBM,adjusted_Pb)

if mode is 'M':
    amplitude=[0,7,100,125,150,15,22,43,49,54,60]

if mode is 'R':
    amplitude=[0,8,100,125,150,17,26,34,43,48,53,58,63,69]


#symmeterize data for fit so we don't get ski-jump :)
#negative_amp=[-value for value in amplitude]
#amplitude_sym=np.concatenate((amplitude,negative_amp))
#Pb_VBM_sym=np.concatenate((Pb_VBM, Pb_VBM))
#Pb_CBM_sym=np.concatenate((Pb_CBM, Pb_CBM))

#lots_o_points = np.linspace(amplitude[0],amplitude[-1],1000, endpoint=True)

#VBMfit = np.polyval(np.polyfit(amplitude_sym, Pb_VBM_sym, 2), amplitude)
#CBMfit = np.polyval(np.polyfit(amplitude_sym, Pb_CBM_sym, 2), amplitude)

#print( "VBM polyfit coefficient for x^2 is: " + str(np.polyfit(amplitude_sym,
#        Pb_VBM_sym, 2)[0]))
#print( "CBM polyfit coefficient for x^2 is: " + str(np.polyfit(amplitude_sym,
#        Pb_CBM_sym, 2)[0]))

np.savetxt('Pb_d_'+mode, np.array((Pb,
           adjusted_Pb)).T, fmt='%.4f')
np.savetxt('VBMCBM_'+mode,np.array((VBM,CBM)).T, fmt='%.4f')
np.savetxt('adjusted_VBMCBM_'+mode,np.array((amplitude,
        VBM,adjusted_CBM)).T,fmt='%.4f')
np.savetxt('Pb_VBMCBM_'+mode, np.array((amplitude,
        Pb_VBM, Pb_CBM)).T, fmt='%.4f')
np.savetxt('all_'+mode, np.array((amplitude, VBM,adjusted_CBM, Pb_VBM, Pb_CBM)).T, fmt='%.4f')
