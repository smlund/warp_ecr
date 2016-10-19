Alfonse N. Pham, Updated Oct. 19, 2016

Initial ion distribution data files serialized with cPickle. The
ions are uniformly distributed in position overlaying the plasma
etching pattern on the injection bias disc. Ion velocities follow
a Maxwell-Boltzmann distribution with energy related to Bohm's
criteria for minimum velocities of acoustic waves in a plasma. 
The data should be formatted such that:

(Species name)_(elecIonEnergy[eV])_np(number of particles).dist

The file extension must be '.dist' 

Work-in-progress
 - git push version1 of a script to generate ions
