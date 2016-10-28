#!/bin/bash
# 
# getReqFilesExtraction.sh
#
# bash script to fetch files required for SuSI extraction simulations
# Written by Alfonse N. Pham, updated Oct. 27, 2016
#
# Note: getReqFilesExtraction.sh must be in folder ./ExtractionSuSI 
# for it to work properly. In the future we can write a more intelligent
# script to place the fetched files into the proper folders regardless
# of script location.
#

# Make 'RequiredFiles' directory if not already present
mkdir -p ../RequiredFiles/initDistExtAper
mkdir -p ../RequiredFiles/SuSIExtraction/Electrodes
cd ../RequiredFiles

# Fetch particle distribution at the extraction aperture for 16Oxygen
cd initDistExtAper
wget -nc https://www.dropbox.com/s/f7ywya8jvdk7pmn/SourceSuSI_v1_000_10.end
cd ..

# Fetch axial magnetic field data for SuSI for BminBecr=0.5,0.73,0.8
cd SuSIExtraction
wget -nc https://www.dropbox.com/s/s8f9nv6qxf4xe3p/ExtracSolenB_BminBecr0.5.csv
wget -nc https://www.dropbox.com/s/05uxg1y3a2kdxek/ExtracSolenB_BminBecr0.8.csv
wget -nc https://www.dropbox.com/s/tap5fafgqebw9kf/ExtracSolenB_BminBecr0.73.csv

# Fetch extraction electrode geometries for diam=8,10,12,14mm
cd Electrodes
mkdir -p SuSI_ExtracAper8mm
cd SuSI_ExtracAper8mm
wget -nc https://www.dropbox.com/s/t9x109efjsl4gv2/Electrodes.dat
wget -nc https://www.dropbox.com/s/bf9susieh9en2r5/Electrodes.wob
cd ..

mkdir -p SuSI_ExtracAper10mm
cd SuSI_ExtracAper10mm
wget -nc https://www.dropbox.com/s/5jxd4j9kpad6n55/Electrodes.dat
wget -nc https://www.dropbox.com/s/huptwkmr2krf3ld/Electrodes.wob
cd ..

mkdir -p SuSI_ExtracAper12mm
cd SuSI_ExtracAper12mm
wget -nc https://www.dropbox.com/s/zxzygdu8qrgqpvi/Electrodes.dat
wget -nc https://www.dropbox.com/s/n4u3zphisywbwd9/Electrodes.wob
cd ..

mkdir -p SuSI_ExtracAper14mm
cd SuSI_ExtracAper14mm
wget -nc https://www.dropbox.com/s/kb2rtdsr4w5wjwp/Electrodes.dat
wget -nc https://www.dropbox.com/s/0xjzuouufiqk0fw/Electrodes.wob
cd ..

cd ../..