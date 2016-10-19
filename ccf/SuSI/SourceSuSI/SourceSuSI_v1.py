# ###------------------------------------------------------------------------------#
# Name:        	SourceSuSI_v1.py
# Purpose:     	Track particle distribution from bias disc to extraction aperture
# Author:      	Alfonse Pham
#				Daniel Winklehner
#
# v1: First version pushed on github Oct. 18, 2016 by A. Pham
# 

import warpoptions

warpoptions.parser.add_argument("-n",
								"--neutralization",
								type=int,
								dest="neutr",
								help="Pass the neutralization of the beam\n"
									"as integer in percent (e.g. -n 80).\n"
									"Default: 100% neutralization (e.g. -n 100)",
								default=100)
warpoptions.parser.add_argument("-a",
								"--autorun",
								action="store_true",
								dest="AutoRun",
								help="Invoke autorun to bypass breakpoints that\n"
									"require user response for diagnostic purposes.",
								default=False)

warpoptions.parser.add_argument("-g",
								"--generate",
								type=str,
								dest="genPath",
								help="Path to existing initial particle\n"
									"distribution file.")

# --- import the necessary scripts to run the simulation properly --- #
# --- WARP imports --- #
from warp import *
from warp.lattice.lattice import *
from warp.particles.extpart import *

# --- Python imports --- #
from numpy import *
import shutil, glob
import sys
import cPickle
import os

# --- Offset of particle positions in vertical direction --- #
# ANP Note: Daniel found that a vertical offset is needed for a symmetric pattern at extraction
# For the field map that I produced with Maxwell, this problem was not there.
#offset = 1.0e-3 # m
offset = 0.0

# --- DW's imports and helper functions--- #
execfile('../../lib/Functions-v3.3.py')
execfile('../../lib/physConst.py')


# --- Set comment lines 2 (main) and 1 (aux) and user's name --- #
top.pline2		= "SuSI Bias Disc to Extraction Aperture"
top.pline1		= "Oxygen, SuSI_B_3D_BminBecr0.80_step1mm.dat"
top.runmaker	= "Alfonse N. Pham"
SimulationPath  = "./results"
# --- Now define the parameters that can be overridden by command line options --- #
autorun			= False		# Stops the execution of the script a certain points and waits for acknowledgement by user.
neut 			= 100		# Beam neutralization (in percent, 100 percent neutralization = virtually no current, hence no space charge effects)

# --- Now override the default values with command line options --- #
if warpoptions.options.AutoRun is not None:
	autorun=warpoptions.options.AutoRun
if warpoptions.options.neutr is not None:
	neut=warpoptions.options.neutr
if warpoptions.options.genPath is not None:
	fname=warpoptions.options.genPath
else:
	fname = ("../RequiredFiles/initDistBiasDisc/16Oxygen_3eV_np100000.dist")
	print ("fname = "+fname)

# --- Change neutralization percentage into a multiplicative factor --- #
neut=(100-float(neut))/100
if abs(neut)<0.0001: neut=1.0e-36 # Note: Warp can't handle 0 current! So set neut to very small value instead.

# --- Invoke setup routine for the plotting (mandatory) --- #
setup()

# --- get root output file name and run number (this is just for convenience) --- #
RunPrefix = arraytostr(top.runid)
RunNumber = setup.pnumb

# --- Create a directory for all the output files ---------------------------- #
SavePath = os.path.join(SimulationPath,RunPrefix+"_"+RunNumber)
#Function mkdir_p() is /lib/Functions-v3.3.py
if not os.path.isdir(SavePath): mkdir_p(SavePath)
# ---------------------------------------------------------------------------- #

# --- Set solver geometry --- #
#w3d.solvergeom=w3d.XYgeom

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# --- Particle/beam parameters initialization, loading and injection --- #

# --- Set number of species --- #
# Note: the function setnspecies also reallocates all the necessary arrays
# So there is no need to call gchange("InPart")...
setnspecies(8)
npSet=100000

# --- Set general flags --- #
top.lrelativ = False                # Relativistic treatment of particles on/off

# --- Load the particle distribution from a .dist file --- #

try:
	FILE = open(fname,'rb')

except IOError:
	print "Couldn't open the particle distribution file. Aborting simulation"
	raise SystemExit

try:
	DISTRIBUTION = cPickle.load(FILE)

except EOFError:
	print("Couldn't unpickle the particle distribution file, maybe not a cPickle file? Aborting simulation")
	raise SystemExit

FILE.close()

# --- Perform a sanity check on the array sizes --- #
# --- Check if array has the right shape --- #
if len( DISTRIBUTION["x"].shape) == 1:
	DISTRIBUTION["x"] = array([DISTRIBUTION["x"]])
	DISTRIBUTION["y"] = array([DISTRIBUTION["y"]])
	DISTRIBUTION["z"] = array([DISTRIBUTION["z"]])
	DISTRIBUTION["vx"] = array([DISTRIBUTION["vx"]])
	DISTRIBUTION["vy"] = array([DISTRIBUTION["vy"]])
	DISTRIBUTION["vz"] = array([DISTRIBUTION["vz"]])
	print "Reshaped the Particle input arrays!"
	if autorun == False:
		breakpoint()

elif len(DISTRIBUTION["x"].shape) > 2 or len(DISTRIBUTION["x"].shape) < 1:
	print "Serious trouble with particle distribution from file . Aborting simulation!"
	raise SystemExit


# --- Set number of species --- #
# nspecies in WARP will be set by determining structure of initial particle distribution
setnspecies(len(DISTRIBUTION['np']))

# --- Check if number of particles per species is correct --- #
if	(DISTRIBUTION["x"].size == DISTRIBUTION["y"].size) and\
	(DISTRIBUTION["y"].size == DISTRIBUTION["z"].size) and\
	(DISTRIBUTION["z"].size == DISTRIBUTION["vx"].size) and\
	(DISTRIBUTION["vx"].size == DISTRIBUTION["vy"].size) and\
	(DISTRIBUTION["vy"].size == DISTRIBUTION["vz"].size) and\
	(DISTRIBUTION["vz"].size == top.ns * npSet):
		pass

else:
	print "Serious trouble with particle distribution from file.\nThe number of particles is probably wrong.\nAborting simulation!"
	raise SystemExit

# --- Set basic beam parameters --- #
#   Note: Only the multispecies parameters (*_s) are set. The single species beam
#   parameters (i.e. without the subscript _s) are set automatically to the values
#   of the first species if not specified otherwise.

#top.aion_s = array([39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948],'d') # Argon-40 species masses [AMU]
#top.zion_s = array([3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0],'d') # charge state of Argon-40 (*1.67e-19 C)
top.aion_s = array([15.995, 15.995, 15.995, 15.995, 15.995, 15.995, 15.995, 15.995],'d') # Oxygen-16 species masses [AMU]
top.zion_s = array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0],'d') # charge state of Oxygen-16 (*1.67e-19 C) 
M_Q = top.aion_s/top.zion_s	# mass-to-charge ratio (for convenience)
top.sp_fract = array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],'d') # species weight
# Note: At least in constant current injection mode, the number of injected simulation particles per species will
# be multiplied with the respective sp_fract value!

# Axial beam velocity (m/s) from Bohm Criteria where Eelec=Eion=3eV
# Work-in-progress: read function to take in electron and ion temperatures and populate variable top.vbeam_s
top.vbeam_s = array([6015.148, 6015.148, 6015.148, 6015.148, 6015.148, 6015.148, 6015.148, 6015.148])	
top.ekin_s = array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])             # kinetic energy                (eV)
non_neutr_ibeam	= array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])        # Non neutralized beam current  (mA)
top.ibeam_s =  neut*non_neutr_ibeam*1.0e-3	# total injected beam current times neutralization factor (A)

top.vthperp_s	= ones(top.ns)          # perpendicular temperature 	(eV) TO BE CHANGED DW'11
top.vthz_s		= ones(top.ns)          # longitudinal temperature 		(eV) TO BE CHANGED DW'11

# --- Set initial beam width and convergence/divergence angle --- #
# Note: Since this is set from a distribution file, x.max and y.max are used as beam size
# and the angles are 0, because the x,y,z-velocities are set from the .dist file.
# Note further: There seems to be a bug in WARP, if the first species has a smaller radius
# than the others, some species get cut off in +x direction (seen in RZ simulation). a0_s, b0_s
# both set to maximum beamsize of all species for now.
top.a0_s        = DISTRIBUTION["x"].max(axis=1).max()*1.0e-3	# beam radius in x-direction (m)
top.b0_s        = DISTRIBUTION["y"].max(axis=1).max()*1.0e-3	# beam radius in y-direction (m)
top.ap0_s       = zeros(top.ns, 'd')							# half-angle in x-direction	(rad)
top.bp0_s       = zeros(top.ns, 'd')							# half-angle in y-direction	(rad)

# --- the z of all particles is set to zbeam at initial creation --- #
wxy.lcommonz	= True
top.zbeam       = 0.09 # Bias disc 9cm from starting face of first injection solenoid

# --- Set some conversion factors for cw beam in XY slice mode --- #
# Note: sw corresponds to the formula on p 32 of the manual with Z_len set to 1.
# DO NOT set this manually for RZ and 3D mode! If continuing from a previous run
# ibeam_s and np_s have been changed accordingly, so use same sp_frac as last run!
top.pgroup.sw	= top.ibeam_s*top.sp_fract/(top.vbeam_s*echarge*top.zion_s*int(DISTRIBUTION["np"].copy()))	# ratio of real particles to code particles
top.pgroup.sm	= top.aion_s*amu														# real mass/species (kg)
top.pgroup.sq	= top.zion_s*echarge													# real charge/species (C)

# --- Set boundary conditions for particles --- #
top.pbound0		= absorb		# boundary at z=0
top.pboundnz	= absorb		# boundary at z=endrun
top.pboundxy	= absorb		# boundary at r=rmax

# --- Call deriveqty() to calculate some of the input values (mandatory)--- #
#	Note: WARP looks for input values that are set to 0.0 and tries to calculate them
#   From other input values (e.g. ekin from vbeam and vice versa) deriveqty also sets
#   the single species variables from multispecies variables and vice versa by calling
#	species().
derivqty()

print "### Particles loaded ###"
# --- Wait for user if mode is not autorun --- #
if autorun == False:
	breakpoint()


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# --- Mesh and fieldsolver --- #

#rPipe = 0.04
rPipe = 0.05	#ANP: found to be larger in drawings

# --- radius of beampipe (particle absorbing wall) --- #
top.prwall		= rPipe

# --- Set boundary conditions for field solve --- #
# Type of boundary condition at plane z=0, z=end and sides (x,y)
# 0 is constant potential (dirichlet), 1 is zero normal derivative (neumann)
# and 2 is periodic.
w3d.bound0		= 1				# boundary at z=0
w3d.boundnz		= 1				# boundary at z=endrun
w3d.boundxy		= 0				# boundary at r=rmax

# --- Set field grid size --- #
w3d.xmmin		= -rPipe*1.
w3d.xmmax		= +rPipe*1.
w3d.ymmin		= -rPipe*1.
w3d.ymmax		= +rPipe*1.
##w3d.zmmin       = -0.25
##w3d.zmmax       = 0.251

# --- set number of meshpoints nx, ny and nz --- #
# Note: For FFT this has to be a power of 2, for multigrid solver
# this has to be an even number!
w3d.nx			= 512
w3d.ny			= 512
##w3d.nz          = 32

# --- Length of one slice in z direction (m) --- #
# Note: As a rule of thumb, this should be smaller than 1.0e-3
wxy.ds          = 5.0e-5

print "### Mesh and Fieldsolver variables set ###"
# --- Wait for user if mode is not autorun --- #
if autorun == False:
	breakpoint()


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# --- Lattice elements --- #
B_Path=("../RequiredFiles/SuSI3DFieldMaps/SuSI_B_3D_BminBecr0.50_step1mm.dat")
#B_Path=("../RequiredFiles/SuSI3DFieldMaps/SuSI_B_3D_BminBecr0.73_step1mm.dat")
#B_Path=("../RequiredFiles/SuSI3DFieldMaps/SuSI_B_3D_BminBecr0.80_step1mm.dat")

with open(B_Path, 'rb') as FILE: SuSI_B = cPickle.load(FILE)

print "zstart of B-field: %.2f mm"%(SuSI_B["zs"]*1000.0)
print "zend of B-field: %.2f mm"%((SuSI_B["zs"] + SuSI_B["nz"] * SuSI_B["dz"])*1000.0)

if not autorun:	breakpoint()

# --- Add the SuSI B field as gridded lattice element in WARP --- #
addnewbgrd(	zs	= SuSI_B["zs"],
			ze	= SuSI_B["zs"] + SuSI_B["nz"] * SuSI_B["dz"],
			xs	= SuSI_B["xs"],
			ys	= SuSI_B["ys"],
			sf	= 0.0,
			sc	= 1.0,
			dx	= SuSI_B["dx"],
			dy	= SuSI_B["dy"],
			bx	= SuSI_B["Bx"],
			by	= SuSI_B["By"],
			bz	= SuSI_B["Bz"],
			nx	= SuSI_B["nx"],
			ny	= SuSI_B["ny"],
			nz	= SuSI_B["nz"])

# --- Delete SuSI_B to free memory --- #
del SuSI_B

print "### Lattice elements initialized ###"
# --- Wait for user if mode is not autorun --- #
if autorun == False:
	breakpoint()

# -- Generate a colorlist to set a specific color for each species --- #
clrs=range(0,200,int(200/top.ns))

# --- Set palette for gist drawing --- #
palette('rainbow.gp')

# --- Define a one-line status line for xy slice mode to keep track of simulation --- #
def statusline(frequency=20):
	if top.it%frequency==0:
		z = top.zbeam
		print "Steps: %i, z = %.4f m, Particles alive = %i"%(top.it, z, top.nplive)

# --- install the statusline function to be called after each step --- #
installafterstep(statusline)
print "### Plotting and data saving set up ###"
# --- Wait for user if mode is not autorun --- #
if autorun == False:
	breakpoint()

# --- Generate the PIC code (allocate storage, load ptcls, t=0 plots,... ) --- #
package("wxy") # select package
generate()     # Generate the PIC code

# --- For each species call WARPs convenience function addparticles --- #
for i in range(top.ns):
	addparticles(x = DISTRIBUTION["x"][i,:].copy(),
		     y = DISTRIBUTION["y"][i,:].copy()+ offset,
		     z = top.zbeam,
		     vx	= DISTRIBUTION["vx"][i,:].copy(),
		     vy	= DISTRIBUTION["vy"][i,:].copy(),
		     vz	= DISTRIBUTION["vz"][i,:].copy(),
		     js	= i,
		     resetmoments = True )

loadrho()
fieldsolve(-1)

# Extraction aperture at z=0.494m referenced to 
ZTargetStep	= int((0.494-0.09)/wxy.ds)
z		= 0.09
nstep		= 0
Envelopes 	= XY_init_envelope(top.ns, ZTargetStep)

filename = os.path.join(SavePath, RunPrefix+"_"+RunNumber)

while z < 0.495:
	z = top.zbeam

	# save initial particle data
	if nstep == 0:
		XY_save_particles(filename = filename+"_00.dat")

	# save envelopes along the way
	if nstep <= ZTargetStep:
		XY_write_envelope(Envelopes, nstep)
	if nstep == int(ZTargetStep*(1.0/10.0)):
		XY_save_particles(filename = filename+"_01.dat"	)
	if nstep == int(ZTargetStep*(2.0/10.0)):
		XY_save_particles(filename = filename+"_02.dat"	)
	if nstep == int(ZTargetStep*(3.0/10.0)):
		XY_save_particles(filename = filename+"_03.dat"	)
	if nstep == int(ZTargetStep*(4.0/10.0)):
		XY_save_particles(filename = filename+"_04.dat"	)
	if nstep == int(ZTargetStep*(5.0/10.0)):
		XY_save_particles(filename = filename+"_05.dat"	)
	if nstep == int(ZTargetStep*(6.0/10.0)):
		XY_save_particles(filename = filename+"_06.dat"	)
	if nstep == int(ZTargetStep*(7.0/10.0)):
		XY_save_particles(filename = filename+"_07.dat"	)
	if nstep == int(ZTargetStep*(8.0/10.0)):
		XY_save_particles(filename = filename+"_08.dat"	)
	if nstep == int(ZTargetStep*(9.0/10.0)):
		XY_save_particles(filename = filename+"_09.dat"	)

	if nstep == ZTargetStep:
		XY_save_particles(filename = filename+"_10.dat")
		XY_save_envelope(Envelopes, filename = filename+".env"	)

	# advance 1 step
	nstep=nstep+1
	step(1)

# Couldn't find a good way of telling warp where to put the cgm files
# so this is the work around.
files = glob.iglob(os.path.join(os.getcwd(),'*.cgm'))
for i in files:
    if os.path.isfile(i):
        shutil.move(i,os.path.join(os.getcwd(),'results'))

files = glob.iglob(os.path.join(os.getcwd(),'*.cgmlog'))
for i in files:
    if os.path.isfile(i):
        shutil.move(i,os.path.join(os.getcwd(),'results'))