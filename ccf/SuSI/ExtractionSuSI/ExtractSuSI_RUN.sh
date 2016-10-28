#!/bin/bash
# Multispecies run for SuSI Extraction Simulation

# Original Simulation invocation script in python
SIM_NAME=ExtractSuSI_v1.py

# Open file containing path to the latest WARP executable
WARP_EXEC=python

RUNID=oxygen_BminBecr0.80

# Move data over to DATAPATH for cataloging
DATAPATH=data000

mkdir -p ./${DATAPATH}

# Change "--current" according to measurements
${WARP_EXEC} ${SIM_NAME} -a --current=603.7031 --solvemode=RZ --pullpos=60 --runid=${RUNID} -a
${WARP_EXEC} ${SIM_NAME} -a --current=603.7031 --solvemode=triangle --pullpos=60 --runid=${RUNID} -a

mv *.cgm ./${DATAPATH}
mv *.cgmlog ./${DATAPATH}
mv ./Results ./${DATAPATH}
mv ./RelaxedPotentials ./${DATAPATH}

exit 0