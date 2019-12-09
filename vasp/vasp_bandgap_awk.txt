#!/usr/bin/awk -f 
# Pass me an 'OUTCAR' and I'll extract fundamental gaps at all k-points.

#Looks like:
# k-point    14 :      -0.3333   -0.3333    0.3333
#  band No.  band energies     occupation 
#      1     -18.3955      2.00000
#      2     -18.3927      2.00000

/k-point/ && (NF==6) {
    kpoint=$2
    kx=$4
    ky=$5 
    kz=$6
    printf("k-point: %d kx: %f ky: %f kz: %f ",kpoint,kx,ky,kz)

    getline # Ignore this line 'band No. band energies occupation'
    while (getline && $3 > 1.0 ) #while more than half an electron
    {
        occ=unocc # juggle temporary variables
        unocc=$2 # band energy
    }

    printf("VBM: %f CBM: %f Bg: %f\n",occ,unocc,unocc-occ)
}
