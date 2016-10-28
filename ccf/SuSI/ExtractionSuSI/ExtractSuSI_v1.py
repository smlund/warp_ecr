# ###------------------------------------------------------------------------------#
# Name:         ExtractSuSI_v1.py
# Purpose:      Using igun methods to resolve the plasma meniscus given electrode
#               potentials and geometries. Multi-species particle distribution are
#               then tracked from the extraction aperture to z=150mm 
# Author:       Alfonse Pham
#               Daniel Winklehner
#
# v1: First version pushed on github Oct. 27, 2016 by A. Pham
# 
# -*- coding: utf-8 -*-
"""
"""
# Command Line Option Parser (must be imported before "warp"):
import warpoptions

# add commandline option to set pullerposition
warpoptions.parser.add_argument("-z",
                                "--pullpos",
                                type=int,
                                dest="PullPos",
                                help="Puller position (integer, mm)",
                                metavar="PULLERPOS")

# add commandline option to set solvemode
warpoptions.parser.add_argument("-s",
                                "--solvemode",
                                type=str,
                                dest="SolveMode",
                                help="input mode (RZ for rz mode, IGUN for IGUN"
                                     " like mode, circle and triangle for 3D"
                                     " modes with transv. temperature)",
                                metavar="SOLVEMODE")

# add commandline option to set runid
warpoptions.parser.add_argument("-r",
                                "--runid",
                                type=str,
                                dest="runID",
                                help="set run ID that will be prefixed to all "
                                     "output files (Max. 40 characters).",
                                metavar="RUNID")

# Flag for automated run
warpoptions.parser.add_argument("-a",
                                "--autorun",
                                action="store_true",
                                dest="AutoRun",
                                help="Turn off waiting for the user at certain "
                                     "points\n as well as some of the "
                                     "unnessecary outputs.",
                                default=False)

# Total beam current
warpoptions.parser.add_argument("--current",
                                type=float,
                                dest="Current",
                                help="Total beam current from source (uA)",
                                default=4800)

# Path of user input particle distribution files with extension .end
warpoptions.parser.add_argument("-i",
                                "--initpath",
                                type=str,
                                dest="initPath",
                                help="Path to existing particle distribution\n"
                                    "file at the extraction aperture.")

# User input puller voltage in units of V
warpoptions.parser.add_argument("--pullpot",
                                type=float,
                                dest="pullerVoltage",
                                help="Puller electrode potential (V).",
                                default=-250.0)

# User input extraction voltage in units of V
warpoptions.parser.add_argument("--extpot",
                                type=float,
                                dest="extVoltage",
                                help="Extraction electrode potential (V).",
                                default=20000.0)

# --- Import Scripts --- #
from warp import *
from warp.field_solvers.generateconductors import *
from warp.particles.particlescraper import *
from warp.run_modes.egun_like import *
from warp.particles.extpart import *
import numpy as np
import sys
import cPickle
from numpy.random import *

# --- Set up some variables for capturing the error and # of iterations --- #
from cStringIO import StringIO

mystdout = StringIO()
old_stdout = sys.stdout
mgiters = 100
mgerror = 1.0

# --- USER INPUT --- #
INPUT = {"mode": "RZ",  # Options: "RZ" or one of three 3D options: "IGUN" = randomly spaced circular beam, no
                        # transversal temp, "circle" = same with boltzmann velocity distr., "triangle" = from files.
         "NumPart1": 10000,     # Number of particles in RZ run
         "NumPart2": 10000,     # Number of particles in 3D run (load_particles overrides this)
         "V_ext": 20000.0,      # extraction voltage (V)
         "V_puller": -250.0,    # puller voltage (V)
         "rpipe": 0.05,         # extent in r-direction (m)--beam pipe radius
         "plasmawidth": 0,      # start of plasma in absolute sim coordinates (between 0 and 7 mm)
         "exppull": 65.0,       # puller setting (distance from plasma aperture, min = 15 mm)
         "pline2": "SuSI source",
         "pline1": "Oxygen BminBecr=0.8",
         "runmaker": "Alfonse Pham",
         "runid": "OxygenRun",
         "particle_path": "../RequiredFiles/initDistExtAper/SourceSuSI_v1_000_10.end",
         "normalize": True,  # normalize current in a brute force way to total extracted beam current?
         #"Itot": 0.000638283, #Oxygem BminBecr=0.5
         #"Itot": 0.000604049, #Oxygem BminBecr=0.8
         "Itot": 0.0006037031, #Oxygem BminBecr=0.8 new see csd spreadsheet
         # total (drain) current in amps, used with normalization ... Higher Current poses problems! DW'10
         "axialtemp": 0.0,  # axial temperature (eV)
         "perptemp": 0.0,  # perpendicular temperature (eV)
         "electrontemp": 3.0,  # screening electron temperature (eV)
         "iontemp": 3.0,  # ion temperature (eV)
         "dx": 1.0e-3,  # Cell size in x (m)
         "dy": 1.0e-3,  # Cell size in y (m)
         "dz": 1.0e-3,  # Cell size in z (m)
         "be_zmax": 0.015,  # max z to include Boltzmann electrons (m)
         "zpend": 0.15,  # position to extract particle information at
         "sourceaperture": 6.0e-3,  # extraction aperture radius (m)
         "sourceradius": 5.9e-3,  # radius of injection source (i.e. 'beam')
         "accuracy": 1.0e-5,  # Accuracy of the field solver (V) (1.0e-5 is WARP default for rz)
         "tri_flag": False,
         # Flag whether the sourceradius should be recalculated to match area of aperture filling triangle
         "larger_flag": False  # Make radius larger again (halfway between full aperture and triangle area)
         }

IONS = {"aion": np.ones(8, 'd') * 15.995,  # Ion species masses (amu)
        "zion": np.linspace(3, 12, 10),  # Ion species charge states (multiples of e)
        #"current": np.array([12.090, 22.969, 38.317, 84.547, 121.68, 226.535, 92.979, 39.168], 'd') * 1.0e-6,  # Ion species currents (A), BminBecr=0.5
        "current": np.array([10.706, 18.793, 26.586, 61.831, 64.118, 190.2, 145.805, 85.664], 'd') * 1.0e-6,  # Ion species currents (A), BminBecr=0.8
        "ns": 8,  # Number of species (int!)
        "np": np.ones(8, 'i') * INPUT["NumPart1"]  # Number of particles (loading from file overrides this)
        }

# --- Set filenames and paths ------------------------------------------------ #
# Set path to simulation directory
SimulationPath = os.getcwd()

# Set path to field files
#FieldPath = "/home/alfonse/research_workspace/warpECRIS/files/ExtractionSuSI"
FieldPath = "../RequiredFiles/SuSIExtraction"


# Set path to initial particle distribution files
if warpoptions.options.initPath is not None:
    ParticlePath=warpoptions.options.initPath
    print ("ParticlePath = "+fname)
else:
    ParticlePath=INPUT['particle_path']

# Set path to library files containing shared functions
LibPath = "../lib"

# Set Filename containing important subroutines
Functions_FN = "Functions-v3.3.py"
Functions_Path = os.path.join(LibPath,Functions_FN)

# Path to source extraction axial field data. Uncomment one.
# We can make this prettier in the future.
#Solenoid_FN = "ExtracSolenB_BminBecr0.5.csv"
#Solenoid_FN = "ExtracSolenB_BminBecr0.73.csv"
Solenoid_FN = "ExtracSolenB_BminBecr0.8.csv"
Solenoid_Path = os.path.join(FieldPath,Solenoid_FN)

# Path to electrode geometries and potentials data files for
# various extraction aperture diamter. Uncomment one.
#ElectrodePath = ".../RequiredFiles/SuSIExtraction/Electrodes/SuSI_ExtracAper8mm"
#ElectrodePath = ".../RequiredFiles/SuSIExtraction/Electrodes/SuSI_ExtracAper10mm"
#ElectrodePath = "../RequiredFiles/SuSIExtraction/Electrodes/SuSI_ExtracAper12mm"
ElectrodePath = (FieldPath+"/Electrodes/SuSI_ExtracAper12mm")
#ElectrodePath = ".../RequiredFiles/SuSIExtraction/Electrodes/SuSI_ExtracAper14mm"
Electrode_FN = "Electrodes.dat"
Electrode_Path = os.path.join(ElectrodePath,Electrode_FN)
# ---------------------------------------------------------------------------- #
# --- END USER INPUT --- #

# --- Initialize the helper functions that are in a different file --- #
execfile(Functions_Path)

# --- Check for necessary output folder and create if necessary --- #
CWD, RESULTS_PATH, RELAXED_PATH = SuSI_EX_Output_Folders(SimulationPath)

# --- Override default runid (name of input file stripped of the '.py') --- #
# (the name of the run and all output files will have the runid as prefix)
if warpoptions.options.runID is not None:
    if len(warpoptions.options.runID) > 40:
        INPUT["runid"] = str(warpoptions.options.runID)[:40]
    else:
        INPUT["runid"] = warpoptions.options.runID

# --- Override Puller setting --- #
if warpoptions.options.PullPos is not None:
    INPUT["exppull"] = float(warpoptions.options.PullPos)

# --- Override Solvemode setting --- #
if warpoptions.options.SolveMode in ["RZ", "IGUN", "circle", "triangle"]:
    INPUT["mode"] = warpoptions.options.SolveMode

# --- Set autorun option --- #
autorun = warpoptions.options.AutoRun

# --- Get current from options --- #
# Note: current in option parser is given in uA
INPUT["Itot"] = warpoptions.options.Current * 1.0e-6

# --- Set comment lines, user's name.
top.pline2 = INPUT["pline2"]
top.pline1 = INPUT["pline1"]
top.runmaker = INPUT["runmaker"]
top.runid = INPUT["runid"] + "_" + INPUT["mode"]

# --- Invoke setup routine for the plotting --- #
setup()

# ---  Get root output file and run number (this is just for convenience)  --- #
RunPrefix = arraytostr(top.runid)
RunNumber = setup.pnumb


# --- Load particle data from file and override IONS if mode == triangle --- #
if INPUT["mode"] in ["RZ", "triangle"]:
    DISTRIBUTION = load_particles(ParticlePath)
    print
    print "Loaded particle distribution from file:"
    print "================================================================================"
    print "Masses:", DISTRIBUTION["M"]
    print "Chargestates:", DISTRIBUTION["Q"]
    print "Number of Particles:", DISTRIBUTION["np"]
    IONS["aion"] = DISTRIBUTION["M"]
    IONS["zion"] = DISTRIBUTION["Q"]
    IONS["ns"] = DISTRIBUTION["ns"]
    IONS["np"] = DISTRIBUTION["np"]

    print
    print "Beam vz (file):"
    for i in range(DISTRIBUTION["ns"]):
        print np.mean(DISTRIBUTION["vz"][i][:DISTRIBUTION["np"][i]])

    NumPart = np.max(DISTRIBUTION["np"])

    if len(IONS["current"]) != DISTRIBUTION["ns"]:
        print
        print "Error: Number of species in loaded distribution doesn't match currents! Aborting"
        exit(1)

# --- Set geometry of the fieldsolver --- #
w3d.solvergeom = w3d.RZgeom

if INPUT["mode"] == "RZ":
    NumPart = INPUT["NumPart1"]
else:
    NumPart = INPUT["NumPart2"]

if INPUT["normalize"]:  # normalize current if desired

    IONS["current"] = IONS["current"] * INPUT["Itot"] / sum(IONS["current"])

# --- Extraction System and Source Setup --- #
runlen = 0.15  # End of the simulation (m)

if INPUT["tri_flag"]:

    INPUT["sourceradius"] = \
        np.sqrt(3 * np.sqrt(3.0) * INPUT["sourceaperture"] ** 2 / 4.0 / np.pi)

    if INPUT["larger_flag"]:
        INPUT["sourceradius"] = \
            (INPUT["sourceradius"] + INPUT["sourceaperture"]) * 0.5  # Test case!

currentdensity = IONS["current"] / (np.pi * INPUT["sourceradius"] ** 2)  # current density (A/m^2)
IONS["vbeam"] = np.sqrt(echarge * (INPUT["electrontemp"] + INPUT["iontemp"]) / (
IONS["aion"] * amu))  # inital beam velocity (acoustic speed of the ions - Bohm criterion) DW'10

#print
#print "Beam vz (Bohm criterion)", IONS["vbeam"]

iondensity = currentdensity / IONS["vbeam"]  # charge/meter
ndens = iondensity / (echarge * IONS["zion"])  # ions/meter
w3d.iondensity = sum(iondensity)  # total ions/meter
w3d.electrontemperature = INPUT["electrontemp"]  # set electron temperature
INPUT["iondensity"] = w3d.iondensity[0]

# --- The expression from the IGUN documentation is used, modified for multiple species.
# --- See Write-up DW'09, phiwall in (V)
phiwall = INPUT["electrontemp"] * (log(sum(ndens * IONS["zion"])) - log(sum(ndens * IONS["zion"] \
    * np.sqrt(2.0 * np.pi * 0.511 / (938.0 * IONS["aion"]) 
    * (1 + INPUT["iontemp"] / INPUT["electrontemp"])))))  # updated by DW'09
w3d.plasmapotential = INPUT["V_ext"] + phiwall
INPUT["Vplasma"] = w3d.plasmapotential[0]

# --- Fix grid cell sizes --- #
# for accurate simulation cell sizes must be <=1.00 mm
w3d.dx = INPUT["dx"]  # Cell size in x [m]
w3d.dy = INPUT["dy"]  # Cell size in y [m]
w3d.dz = INPUT["dz"]  # Cell size in z [m]

# -- range over which Boltzmann electrons are included --- #
w3d.xbemin = 0.0  # min radius [m]
w3d.xbemax = INPUT["rpipe"]  # max radius [m]
w3d.zbemin = INPUT["plasmawidth"]  # min z [m]
w3d.zbemax = INPUT["be_zmax"]  # max z [m]

# --- Set basic beam parameters --- #
top.a0 = INPUT["sourceradius"]  # initial beam width in x
top.b0 = INPUT["sourceradius"]  # initial beam width in y
top.ap0 = 0.0  # initial beam envelope vx/vz ?
top.bp0 = 0.0  # initial beam envelope vy/vz ?
top.emit = 0.0  # initial perp emittance (rms-edge)
top.ekin = 0.0  # beam kinetic energy (if 0 it will be calculated by derivqty() for every species later)
top.lrelativ = False  # flag for relativistic effects

# --- Multiple Species --- #
setnspecies(IONS["ns"])  # set number of species
top.aion_s = IONS["aion"]  # ion masses [amu]
top.zion_s = IONS["zion"]  # ion charges [*1.67e-19 C]
top.vbeam_s = IONS["vbeam"]  # axial beam velocity [m/s]
top.ibeam_s = IONS["current"]  # total injected beam current

#print
#print "Currents:", top.ibeam_s, "A"
#print "Total current:", sum(top.ibeam_s), "A"
#breakpoint()

# --- Distinguish between RZ and 3D run for particle injection --- #
if INPUT["mode"] == "RZ":

    #particles = createRZbeam(NumPart=min(DISTRIBUTION["np"]),
    particles = createRZbeam(NumPart=NumPart,
                             ns=DISTRIBUTION["ns"],
                             radius=INPUT["sourceradius"])

    installuserparticlesinjection(rzInjection)

else:

    installuserparticlesinjection(myinjection)

w3d.l_inj_user_particles = True  # User defines x,y,z of injected particles
w3d.l_inj_user_particles_v = True  # User defines vx,vy,vz of injected particles

##if INPUT["mode"] == "RZ":
##
##	top.vthz_s = np.sqrt(INPUT["axialtemp"]*2.0*echarge*top.zion_s/top.aion_s/amu)	# axial thermal spread by species [m/s]
##	top.vthperp_s = np.sqrt(INPUT["perptemp"]*echarge*top.zion_s/top.aion_s/amu) # perpendicular thermal spread by species [m/s]
##
##else:
if 1:
    # --- Because of the way the injection is done, the particle velocities
    # --- are multiplied by these quantities so they must be set to a nonzero
    # --- value if a finite velocity spread is desired.
    top.vthperp = 1.0
    top.vthz = 1.0

    for i in range(top.ns):
        top.vthperp_s[i] = 1.0
        top.vthz_s[i] = 1.0

M_Q = top.aion_s / top.zion_s  # mass-to-charge ratio
derivqty()  # calculate ekin from vbeam or vice versa

w3d.l4symtry = True  # Set 4-fold symmetry
w3d.l2symtry = False  # No 2-fold symmetry

# --- Set boundary conditions for field solve --- #
w3d.bound0 = dirichlet  # boundary at z=0
w3d.boundnz = neumann  # boundary at z=endrun
w3d.boundxy = neumann  # boundary at r=rmax

# --- Set boundary conditions for particles --- #
top.pbound0 = absorb  # boundary at z=0
top.pboundnz = absorb  # boundary at z=endrun
top.pboundxy = absorb  # boundary at r=rmax (Newly added)
top.prwall = INPUT["rpipe"]  # Radius of cylindrical wall that absorbs particles

# --- Set field grid size --- #
w3d.xmmin = -INPUT["rpipe"] * 1.0  # Mesh lower limit in x
w3d.xmmax = +INPUT["rpipe"] * 1.0  # Mesh upper limit in x
w3d.ymmin = -INPUT["rpipe"] * 1.0  # Mesh lower limit in y
w3d.ymmax = +INPUT["rpipe"] * 1.0  # Mesh upper limit in y
w3d.zmmin = INPUT["plasmawidth"]  # Mesh lower limit in z
w3d.zmmax = runlen  # Mesh upper limit in z

# --- Set grid borders if we have symmetry on --- #
if w3d.l4symtry: w3d.xmmin = 0.0
if w3d.l2symtry or w3d.l4symtry: w3d.ymmin = 0.0

# --- Reset nx, ny, nz for FFT --- #
w3d.nx = goodnn(nint((w3d.xmmax - w3d.xmmin) / w3d.dx))
w3d.dx = (w3d.xmmax - w3d.xmmin) / (w3d.nx * 1.0)
w3d.xmmax = w3d.xmmin + w3d.nx * w3d.dx
w3d.ny = goodnn(nint((w3d.ymmax - w3d.ymmin) / w3d.dy))
w3d.dy = (w3d.ymmax - w3d.ymmin) / (w3d.ny * 1.0)
w3d.ymmax = w3d.ymmin + w3d.ny * w3d.dy
w3d.nz = goodnn(nint((w3d.zmmax - w3d.zmmin) / w3d.dz))
w3d.dz = (w3d.zmmax - w3d.zmmin) / (w3d.nz * 1.0)
w3d.zmmax = w3d.zmmin + w3d.nz * w3d.dz

# --- Timestep size and scaling --- #
ekininit = top.aion_s[1] * amu * top.vbeam_s[1] ** 2 / 2 / echarge  # initial beam kinetic energy of species with min(M_Q)
vmax = np.sqrt(2.0 * echarge * (INPUT["V_ext"] + ekininit - INPUT["V_puller"]) / (amu * min(M_Q)))  # ~max beam velocity
top.dt = 0.3 * w3d.dz / vmax  # set dt to .3*dz/(max beam velocity)
top.pgroup.dtscale[0:top.ns] = np.sqrt(M_Q) / np.sqrt(min(M_Q))  # scale times steps for all other species

# --- Specify injection of the particles --- #
top.npmax = 0  # start with 0 particles since all will be injected
top.inject = 1  # injection type: 1)Const Current, 2)Child-Langmuir s.c. limited, 3) Gauss s.c. limited
top.injctspc = 1000000  # Extra particle storage space
top.rinject = largepos  # Source radius of curvature
top.ainject = INPUT["sourceradius"]  # initial beam width in x
top.binject = INPUT["sourceradius"]  # initial beam width in y
top.npinject = NumPart  # Number of particles injected each step Np
top.linj_eperp = False  # See top.v about this
top.linj_enormcl = False
top.linj_efromgrid = True  # See top.v about this
top.vinject = w3d.plasmapotential  # Voltage on the injection source
w3d.l_inj_addtempz_abs = True  # longitudinal thermal velocity is positive
top.zinject = INPUT["plasmawidth"]  # location of injection
top.lvinject = True

# --- SET UP CALCULATION OF MOMENTS --- #
top.lsavelostpart = False  # Line to save lost particle information (false = no saving)
top.lepsaveonce = True  # Line to save particle information only once
# When true, each particle is saved at most only once in each window
# This only works if the windows do not overlap each other!
# --- For now: turn off all calculation of moments --- #
top.ifzmmnt = 0  # turn off calculation of beam frame moments...
top.iflabwn = 0  # turn off calculation of lab window moments...
top.lspeciesmoments = True  # Set if moments should be caclulated for each species (true) or only for the whole beam (false)

clrs = range(0, 200, int(200 / top.ns))  # set a specific color for each species

# --- Set up fieldsolver --- #
# Note: 7 means the multigrid solver
top.fstype = 7
f3d.mgparam = 1.5
f3d.downpasses = 2
f3d.uppasses = 2

# --- Set source solenoid field --- #
load_SuSI_solenoid(Solenoid_Path, autorun=autorun)
derivqty()  # calculate ekin from vbeam or vice versa
# --- Generate the PIC code (allocate storage, load ptcls, t=0 plots, etc.)
package("w3d")  # select package
generate()  # Generate the PIC code

# --- Set accuracy of fieldsolver --- #
# Note: absolute error - greatest difference (in V) of actual field and last field
# Default in WARP is 1.0e-5
frz.mgridrz_accuracy = INPUT["accuracy"]

# --- Set grid to non-moving -> conductors are calculated only once, should save some time
f3d.gridmode = 1

# --- Open plot window
winon()
palette('rainbow.gp')  # set palette for printing

# --- Set up Extraction Electrodes -------------------------------------------------- #
print "Setting up Plasma voltage of %.4e V at z = %.4e" % (w3d.plasmapotential, 1.0e-3 * INPUT["plasmawidth"])
plasma = Box(2 * w3d.xmmax, 2 * w3d.ymmax, 0., w3d.plasmapotential, zcent=1.0e-3 * INPUT["plasmawidth"])
if not autorun: breakpoint()

# Load the conductor meta file
with open(Electrode_Path, "rb") as FILE: electrode_data = cPickle.load(FILE)

voltages = electrode_data["voltages"]
cond_ids = electrode_data["cond_ids"]
names = electrode_data["cond_names"]

if warpoptions.options.pullerVoltage is not None:
    voltages[1] = warpoptions.options.pullerVoltage
else:
    voltages[1] = INPUT["V_puller"]

if warpoptions.options.extVoltage is not None:
    voltages[0] = warpoptions.options.extVoltage
else:
    voltages[0] = INPUT["V_ext"]

# Shift puller and ground according to setting
# There are three conductors in the file, the first (index 0) is the source and
# is not shifted, the other two are puller and ground and are both shifted
zshift = (INPUT["exppull"] - 15.0) * 1e-3
zshifts = [0.0, zshift, zshift]

# Load geometry from *wob file
conds = SRFRVLAfromfile(os.path.splitext(Electrode_Path)[0],
                        voltages, cond_ids, zshifts=zshifts, install=0)

# Select the conductor data from class
conductors = conds.conds[0]

for i in range(len(conds.conds) - 1):
    conductors += conds.conds[i + 1]

# Update plot
fma()
conductors.draw()
refresh()

# Install the conductors and particle scraper
installconductors(plasma)
installconductors(conductors)
scraper = ParticleScraper(conductors=conductors)
zsrcmin = [0.007]
print "Done setting up conductors"
if not autorun: breakpoint()
# ---------------------------------------------------------------------------- #

# --- Get initial fields ----------------------------------------------------- #
f3d.mgmaxiters = 100

if INPUT["mode"] == "RZ":
    # If in RZ mode, calculate fields
    nl = frz.basegrid.nlevels
    frz.basegrid.nlevels = 0
    fieldsolve(-1)
    frz.basegrid.nlevels = nl
    fieldsolve(-1)
    fieldsolve(-1)

else:
    # in other modes, read in the field from file from previous RZ run
    f = PR.PR("./RelaxedPotentials/" + INPUT["runid"] + "_LastRun.wrp")  # should be same as below!
    frz.basegrid.erp = f.ER
    frz.basegrid.ezp = f.EZ
    frz.basegrid.phi[...] = f.pot
    frz.basegrid.rho[...] = f.rho
    f.close()
# ---------------------------------------------------------------------------- #

# --- Install some functions ------------------------------------------------- #
installbeforefs(redirect_stdout)  # install before the fieldsolver
installbeforefs(setplasmarho)  # install before the fieldsolver
installafterfs(get_vcycles)  # install directly after the fieldsolver
installafterstep(plot_xy)
# ---------------------------------------------------------------------------- #
print
print "Currents:", top.ibeam_s, "A"
print "Total current:", sum(top.ibeam_s), "A"
print "ns:", top.pgroup.ns
print "nps", top.pgroup.nps
print "M:", top.pgroup.sm
print "Q:", top.pgroup.sq
print "SW:", top.pgroup.sw
breakpoint()


# ------------------------ START THE SIMULATION RUN -------------------------- #
"""
# First Run - RZ Mode - Solve for fields
# First run iterations to solve symmetric simulation
# Number of iterations at each level SetByUser
# lvariabletimestep=1,fvariabletimestep=0.3 ???
"""
if INPUT["mode"] == "RZ":

    for i in range(1):
        gun(1, ipsave=200000, rhoparam=0.0, maxtime=5.0e-6)
        plotfields(1)

    for i in range(10):
        gun(1, ipsave=200000, rhoparam=0.0, maxtime=5.0e-6,
            current=IONS['current'], currentiz=slice(-100, -1))
        plotfields(1)

    for i in range(5):
        gun(1, ipsave=200000, rhoparam=0.5, maxtime=5.0e-6,
            current=IONS['current'], currentiz=slice(-100, -1))
        plotfields(1)

    ##    gun(1, ipsave=200000, maxtime = 5.0e-6)

    ##    counter = 0
    ##
    ##    while mgiters > 6 or mgerror > INPUT["accuracy"]:
    ##
    ##        if counter >= 10:
    ##
    ##            print "System didn't converge after 10 gun-runs at rhoparam = 0.0!"
    ##            answer = raw_input("Do you want to use the field anyway? (y/n)")
    ##            if answer == 'y': break
    ##            else: exit()
    ##
    ##        gun(1, ipsave=200000, maxtime = 5.0e-6)
    ##        plotfields(1)
    ##
    ##        print "pgroup.npmax: ", top.pgroup.npmax
    ##        counter +=1

    # --- Next write out solved field to file
    f = PW.PW("./RelaxedPotentials/" + INPUT["runid"] + "_LastRun.wrp")
    f.ER = frz.basegrid.erp
    f.EZ = frz.basegrid.ezp
    f.pot = frz.basegrid.phi
    f.rho = frz.basegrid.rho
    f.rz_inf = np.array([w3d.xmmin, w3d.zmmin, w3d.xmmax, w3d.zmmax,
                         w3d.dx, w3d.dz, w3d.nx, w3d.nz])
    f.IONS = IONS
    f.INPUT = INPUT
    f.close()

# --- In RZ Run one more time to collect particle information at z = 0.15 m
# --- In 3D this is the simulation...
zcp = ZCrossingParticles(laccumulate=True, zz=INPUT["zpend"])
gun(1, ipsave=200000, maxtime=5.0e-6)
plotfields(1)
# ---------------------------------------------------------------------------- #

# ------------ Call some functions to collect and save the data -------------- #
# Envelopes
env_fn = INPUT["runid"] + '_' + INPUT["mode"] + '.env'
save_envelope_3D(zstart=INPUT["plasmawidth"],
                 zend=INPUT["zpend"],
                 filename=os.path.join(RESULTS_PATH, env_fn))
print "Final envelope plots saved!"

# Particle data at zcp
end_fn = INPUT["runid"] + '_' + INPUT["mode"] + '.end'
save_particles_3D(INPUT,
                  IONS,
                  filename=os.path.join(RESULTS_PATH, end_fn),
                  zcp=zcp)
print "Final particle data saved!"

# Write textfile with all the information of the run
inf_fn = INPUT["runid"] + '_' + INPUT["mode"] + '.inf'
save_summary_3D(INPUT,
                IONS,
                filename=os.path.join(RESULTS_PATH, inf_fn))
print "Final summary saved!"

# Write file, that contains the filenames for set loading in postprocessor
set_fn = INPUT["runid"] + '_' + INPUT["mode"] + '.set'
write_set_file(filename=os.path.join(RESULTS_PATH, set_fn))
print " Set file saved"
# -----------------------------------------------------------------------------#