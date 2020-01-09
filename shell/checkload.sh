#!/bin/bash

for i in curie bragg planck kelvin
do
echo " "
ssh $i 'hostname'
ssh $i 'uptime'
cores=$(ssh $i 'grep processor /proc/cpuinfo | wc -l')
echo " There are "$cores "cores"
echo " "
done
