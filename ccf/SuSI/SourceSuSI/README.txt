Alfonse N. Pham, Updated Oct. 19, 2016

Once all field and particle distribution files are in their respective
folders, source simulation is executed in command line in various ways:

# Quick run, default initial distribution file used, and breakpoints are ON
python SourceSuSI_v1.py

# Autorun (no breakpoints) with '-a' option, default initial distribution 
# file used
python SourceSuSI_v1.py -a

# Autorun (no breakpoints) with '-a' option, initial distribution file in the
# path expressed in a string after the '-g' option 
python SourceSuSI_v1.py -a -g "../RequiredFiles/SuSI3DFieldMaps"

