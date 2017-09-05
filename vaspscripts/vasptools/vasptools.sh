#!/usr/bin/bash
# set up environment for vasptools
# Germain Vallverdu <germain.vallverdu@univ-pau.fr>
# 15/05/2012

# checkout directory
INSDIR=/home/gvallverdu/Programme

# add scripts directory to PATH
export PATH=$INSDIR/vasptools/scripts/:$PATH

# add modules directory to PYTHONPATH
export PYTHONPATH=$INSDIR/vasptools:$PYTHONPATH

