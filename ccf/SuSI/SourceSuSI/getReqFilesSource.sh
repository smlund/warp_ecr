#!/bin/bash
# 
# getReqFilesSource.sh
#
# bash script to fetch files required for SuSI source simulations
# Written by Alfonse N. Pham, updated Oct. 20, 2016
#
# Note: getReqFilesSource.sh must be in folder ./SourceSuSI for it
# to work properly. In the future we can write a more intelligent
# script place the fetched files into the proper folders regardless
# of script location.
#

# Make 'RequiredFiles' directory if not already present
mkdir -p ../RequiredFiles/SuSI3DFieldMaps
cd ../RequiredFiles

# Fetch initial particle distribution files for 16Oxygen
cd initDistBiasDisc
wget -nc https://www.dropbox.com/s/f9dek6wwwexwsbc/16Oxygen_3eV_np10000.dist
wget -nc https://www.dropbox.com/s/g9llfdgoq42roes/16Oxygen_3eV_np50000.dist
wget -nc https://www.dropbox.com/s/2ggbz8aqazhhlwe/16Oxygen_3eV_np100000.dist
cd ..

# Fetch 3D magnetic field maps for SuSI for BminBecr=0.5,0.73,0.8
cd SuSI3DFieldMaps
wget -nc https://www.dropbox.com/s/74nwn9r9p3sulxl/SuSI_B_3D_BminBecr0.50_step1mm.dat
wget -nc https://www.dropbox.com/s/ihwg0fzkgfoxkxw/SuSI_B_3D_BminBecr0.73_step1mm.dat
wget -nc https://www.dropbox.com/s/cjop14ncqispwsc/SuSI_B_3D_BminBecr0.80_step1mm.dat
cd ..