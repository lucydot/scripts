import argparse
import pandas as pd
from io import StringIO
import csv

def make_POSCAR(iteration,outcar="OUTCAR",headings_file="POSCAR",fileout="POSCAR_"):
    
    # read in data from outcar and extract list of atomic positions at a given ionic step
    data = open(outcar,'r').read()
    data = data.split("Iteration   "+str(iteration),1)[-1].split("total drift,1")[-1].split("total drift")[0].split("Angst)")[-1]
    df = pd.read_csv(StringIO(data),delimiter="\n",names="r",skiprows=2)
    positions = [df.r[i].split()[:3] for i in range(df.r.count()-1)]
    
    # read the POSCAR heading information and write to new fileout
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
        for line in positions:
            fileout.write(' '.join(n for n in line)+'\n')
         
if __name__== '__main__':
    parser = argparse.ArgumentParser(description='extract positions from OUTCAR at particular iteration and create a POSCAR from them')
    parser.add_argument('-i','--iteration', type=int, help = 'Iteration at which you would like to extract atomic positions')
    parser.add_argument('-fi','--filein', default='OUTCAR',type=str, help='name of file (usually OUTCAR) you would like to extract positions from')
    parser.add_argument('-hi','--heading', default='POSCAR',type=str, help='name of file (usually POSCAR) you would like to extract headings from')
    parser.add_argument('-fo','--fileout',default="POSCAR_",type=str, help='name of file you would like to write to. Iteration number will be suffixed')
    args = parser.parse_args()
    print (args.iteration, args.filein, args.heading, args.fileout)
    make_POSCAR(args.iteration,args.filein,args.heading,args.fileout)

    

