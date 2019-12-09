#!/bin/bash -l

for NUM in {-14,-13,-12,-11,-10,-09,-08,-07,-06,-05,-04,-03,-02,-01,000,001,002,003,004,005,006,007,008,009,010,011,012,013,014}
do
   if [ -d DISP_$NUM ]
   then
     echo DISP_$NUM
     path=$path" DISP_$NUM"
   fi
done

echo $path

get_eig.py -p $path -v ../00_d0/DISP_-14 -e 1.5 5.5 
