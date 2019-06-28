#! /bin/bash
 
# Script to extract the atomic positions (xyz) from 
# POSCAR & OUTCAR files of VASP and write them to a 
# file, to be used with Jmol, MolDen and such other
# third party programs
 
# First written: Gowtham, Tue, 22 Dec 2009 11:11:12 -0500
# Last modified: Gowtham, Tue, 22 Dec 2009 16:29:36 -0500
 
# Essential VASP Files
# If missing, quit the program
export POSCAR="POSCAR"
export OUTCAR="OUTCAR"
 
if [[ -e $POSCAR || -e $OUTCAR ]];
then
  echo
  echo "  $POSCAR & $OUTCAR found"
else
  echo
  echo "  ERROR : POSCAR &/or OUTCAR missing"
  echo
  exit
fi
 
# Output Filename 
# If missing, quit the program
export OUTPUT=$1
 
if [ $# != 1 ]
then
  echo
  echo "  ERROR : OUTPUT_FILENAME Missing"
  echo "  Usage : vasp2xyz.sh OUTPUT_FILENAME"
  echo
  exit
fi
 
# Delete scratch files
rm -f Final.tmp.*
rm -f $OUTPUT.xyz
 
# Atoms related information
# Extracted from POSCAR & OUTCAR
 
# Total number of atom types
# Pat Krogel: Tue, 22 Dec 2009 13:13:48 -050
declare -a ATOM_TYPE=(`sed -n '6 p' $POSCAR`)
export ATOM_TYPE_TOTAL=${#ATOM_TYPE[@]}
 
# Total number of atoms
export ATOMS_TOTAL=0
export i=0
while [[ "$i" -lt "$ATOM_TYPE_TOTAL" ]]
do
  export iATOM_TYPE=${ATOM_TYPE[$i]}
  export ATOMS_TOTAL=`expr $ATOMS_TOTAL + $iATOM_TYPE`
  export i=`expr $i + 1`
done
 
# Atom labels
declare -a ATOM_LABEL=(`grep "POTCAR:" $OUTCAR | head -$ATOM_TYPE_TOTAL | awk '{print $3}' | sed -e 's/_.*//g'`)
export ATOM_LABEL_TOTAL=${#ATOM_LABEL[@]}
 
# Print2Screen
if [[ "$ATOM_TYPE_TOTAL" != "$ATOM_LABEL_TOTAL" ]]
then
  echo
  echo "  ERROR : Atom Type Total ($ATOM_TYPE_TOTAL) does NOT"
  echo "          match Atom Label Count ($ATOM_LABEL_COUNT)"
  echo "  Please make sure OUTCAR & POSCAR are from the same"
  echo "  system and calculation"
  echo 
  exit
else
  echo
  echo "  Number of Atom Types   :" $ATOM_TYPE_TOTAL
  echo "  Atom Type              :" ${ATOM_LABEL[@]}
  echo "  Atom Type Count        :" ${ATOM_TYPE[@]}
  echo "  Total Number of Atoms  :" $ATOMS_TOTAL
 
  # ATOM_LABEL_2 is an array of atom labels
  # with each atom label appearing an appropriate
  # number of times, to a total of $ATOMS_TOTAL
  declare -a ATOM_LABEL_2
  export REPEAT=0
  export j=0
 
  while [[ "$j" -lt "$ATOM_TYPE_TOTAL" ]]
  do
    export k=1
 
    while [[ "$k" -le "${ATOM_TYPE[$j]}" ]]
    do
      export REPEAT=`expr $REPEAT + 1`
      export k=`expr $k + 1`
      ATOM_LABEL_2[$REPEAT]=${ATOM_LABEL[$j]}
    done
 
    export j=`expr $j + 1`
  done
# echo ${ATOM_LABEL_2[@]}
# echo ${#ATOM_LABEL_2[@]}
 
  # Record the line numbers where coordinates are written
  # in OUTCAR. Actual coordinates begin 2 lines below
  # the recorded numbers
  grep -n "POSITION" $OUTCAR | sed -e 's/:/ /g' | sed -e 's/POSITION.*//g' > LINE_NUMBERS.tmp
 
  export FRAME_NUMBER=1
  while read START
  do
    export START_HERE=`expr $START + 2`
    export END_HERE=`expr $START_HERE + $ATOMS_TOTAL - 1`
 
    echo "$ATOMS_TOTAL" >> $OUTPUT.xyz
    echo "# Frame Number: $FRAME_NUMBER" >> $OUTPUT.xyz
    sed -n "$START_HERE,$END_HERE p" $OUTCAR | awk '{printf "%12.8f %12.8f %12.8f\n", $1, $2, $3}' > Final.tmp.$FRAME_NUMBER
 
    export l=1
    while read X Y Z
    do
      echo ${ATOM_LABEL_2[$l]} $X $Y $Z > Final.tmp.$FRAME_NUMBER.$l
      awk '{ printf "%-3s %12.8f %12.8f %12.8f\n", $1, $2, $3, $4 }' Final.tmp.$FRAME_NUMBER.$l >> $OUTPUT.xyz
      export l=`expr $l + 1`
    done<Final.tmp.$FRAME_NUMBER
 
    export FRAME_NUMBER=`expr $FRAME_NUMBER + 1`
  done<LINE_NUMBERS.tmp
 
  export FRAMES_TOTAL=`expr $FRAME_NUMBER - 1`
  echo "  Total Number of Frames :" $FRAMES_TOTAL
  echo
fi
 
# Delete scratch files
rm -f LINE_NUMBERS.tmp 
rm -f Final.tmp.*
