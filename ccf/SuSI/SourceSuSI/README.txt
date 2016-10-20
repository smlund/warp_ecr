Alfonse N. Pham, Updated Oct. 19, 2016

Before simulation can be successfully executed, field maps and particle
distribution files must be fetched and placed in their proper path. To do
so, execute bash script 'getReqFilesSource.sh' by typing into the command
line:

./getReqFilesSource.sh

Once all field map and particle distribution files are in their respective
folders, source simulation is executed in command line given in examples:

# Quickrun with breakpoints ON. Default initial distribution 
# file '16Oxygen_3eV_np10000.dist' and default 3D field map file
# 'SuSI_B_3D_BminBecr0.73_step1mm.dat' are used.
python SourceSuSI_v1.py

# Autorun (no breakpoints) with '-a' option. Default initial distribution 
# file '16Oxygen_3eV_np10000.dist' and default 3D field map file
# 'SuSI_B_3D_BminBecr0.73_step1mm.dat' are used.
python SourceSuSI_v1.py -a

# Autorun (no breakpoints) with '-a' option. Non-default initial distribution 
# file in the path expressed in string after '-i' option and default 3D field 
# map file 'SuSI_B_3D_BminBecr0.73_step1mm.dat' is used.
python SourceSuSI_v1.py -a -i "../RequiredFiles/initDistBiasDisc/16Oxygen_3eV_np10000.dist"

# Autorun (no breakpoints) with '-a' option. Non-default 3D field map file in 
# the path expressed in string after '-b' option and default initial distribution 
# file '16Oxygen_3eV_np10000.dist' are used.
python SourceSuSI_v1.py -a -b "../RequiredFiles/SuSI3DFieldMaps/SuSI_B_3D_BminBecr0.73_step1mm.dat"