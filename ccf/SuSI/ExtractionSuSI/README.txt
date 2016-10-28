Alfonse N. Pham, Updated Oct. 27, 2016

Before the extraction simulation can be successfully executed, axial field
and electrode geometry files must be fetched and placed in their proper 
path. To do so, execute bash script 'getReqFilesExtraction.sh' by typing 
into the command line:

./getReqFilesExtraction.sh

For a default run of the WARP extraction codes, run the bash script
'EtractSuSI_RUN.sh' by typing in the command line in the same path as
this README.txt:

./ExtractSuSI_RUN.sh

Read the bash script carefully as it give you the sequence into which the
extraction code should be executed:
- Initially in RZ mode to establish the plasma sheath and store it in the
  folder '${DATAPATH}/RelaxedPotentials'
- Then in Triangle mode, that takes the particle distribution at the
  extraction electrode and track it through the potential established in the
  RZ mode.

Work-in-progress
================
- optimize code so that the igun routine convergence prompt exit of while
  loop. Preliminary code were written to do this, but will need more testing
  to ensure stability in the subsequent version. For now, convergence is
  established by running the code via trial and error.

 
