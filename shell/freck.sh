#!/bin/bash  
# Force Bash, so we have 'let' for addition 

# Freck - "To move swiftly or nimbly."
# http://matadornetwork.com/abroad/20-obsolete-english-words-that-should-make-a-comeback/

VERSION="0.03"
# 0.01 Works!
# 0.02 Uses find + moves [[:upper:]] files to subdir.
# 0.03 If we appears to be in a git repo; checkin breadcrumb [and thus folder]

echo "If you knew Time as well as I do, said the Hatter, you wouldn't talk about wasting IT. It's HIM."

# Following lines keep a running count in a hidden file .count in the present directory...
count=` cat .count `
let count=count+1
echo $count > .count
# Format this as 0002 etc. ; pad to 4 digits with awk/printf
prefix=` echo $count | awk '{printf("%04d",$1)}' `

# Request name for this set of calcs / data
echo "OK; going to bundle & rename to (please enter...):"
echo -n "${prefix}-"
read suffix

# Glob together for use from this point on
newdir="${prefix}-${suffix}"
# Should look like "0004-MY_wicked_calculations"

echo -n "Making directory ${newdir} ..."
mkdir "${newdir}"
echo "Moving files to ${newdir}..."

#mv *.* "${newdir}" # Standard files
find  ./ -maxdepth 1 -name "*.*" -not -name ".*" -type f -print0 | xargs -0 -I {} mv {} "${newdir}"
find  ./ -maxdepth 1 -name "[[:upper:]]*" -type f -print0 | xargs -0 -I {} mv {} "${newdir}"

# find  ./ -maxdepth 1 -name "[[:upper:]]*" -type f   # something like this + then mv
#mv [A-Z] "${newdir}" # VASP and other Fortran programmes ALLCAPS output files (no extension)

# Automatic logging - would like to add a lot more...

echo "Leaving breadcrumbs..."
crumb="${newdir}/`date +%Y-%m-%d_%H%S`_${USER}_at_${HOSTNAME}"
# Filename looks like: 2015-09-11_0002_jarvist_at_Asriel

touch "${crumb}"

cat > ${crumb} << EOF
#Archived with Freck Version ${VERSION}... 
# Still not able to extract history :(
pwd: ` pwd `
date: ` date `
hostname: ` hostname `
EOF

# Doesn't work currently - we're in a sub shell (~:
# Nice things to do with permanently changing .bashrc to adjust history defaults, which might be an idea anyway...
#   HISTFILE="${newdir}"/`date +%Y-%m-%d`_${USER}_at_${HOSTNAME}.dat
#   history -a # Dump history saved in memory (i.e. this shell)


# GIT integration

if [ -d ".git" ]
then
    echo "This appears to be a git repository... checking in breadcrumb."
    git add "${crumb}"
    git commit -m "${newdir} freck auto commit."

    if [ -e ${newdir}/POSCAR ]
    then
        echo "Detected VASP: generating POTCAR.TITEL, commiting along with POSCAR, INCAR, KPOINTS"
        grep "TITEL"  ${newdir}/POTCAR > ${newdir}/POTCAR.TITEL
        git add ${newdir}/POTCAR.TITEL 
        git add ${newdir}/POSCAR 
        git add ${newdir}/INCAR 
        git add ${newdir}/KPOINTS
        git add ${newdir}/CONTCAR
        git commit -m "${newdir} freck VASP autosave inputs"
    fi

fi
