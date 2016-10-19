Alfonse N. Pham, 2016

3D B-Field Map were generated in Ansys Maxwell and converted
to cPickle serialized object for use with the WARP code in
its current format. The serialization improves the read time
of larger data files.

Work-in-progress
 - Write bash script to git pull SuSI 3D Bfield Maps
 - upload Maxwell files to produce arbitrary 3D magnetic 
   field maps with instructions
 - git push a script to convert the Maxwell data format 
   into cPickle format recognized by the WARP simulations.