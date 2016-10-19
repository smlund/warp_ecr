###------------------------------------------------------------------------------#
# Name:        Functions.py
# Purpose:     Supporting functions and submodules for WARP
# Author:      D Winklehner
#              A Pham
#
# Created:     18/04/2011
# License:     no license
#
# v2.1: Added function to read in Damons file for the VENUS source fields and
# import numpy as np from now on (instead of *)
# v3.0: Added functions from VENUS extraction simulation
# v3.1: Added function to make output folders in Extraction sim
#------------------------------------------------------------------------------#

# --- Imports --- #
import termios
import fcntl
import sys
import os
import numpy
import numpy as np
from numpy.random import random
from struct import unpack
# -----------------------------------------------------------------------------#

# ----- Quick Tools ----- #
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

# -----------------------------------------------------------------------------#
# ---------------------------- OUTPUT FUNCTIONS -------------------------------#
# -----------------------------------------------------------------------------#
def XY_save_particles(mode = 0, filename = None, ibeam_old = None, np_old = None):
    """
    XY_save_particles is a function to conveniently save particles at the
    current XY slice into a cPickled file. It can later be read out with the
    Postprocessor
    Note: Works only in Slice Code simulations.
    mode = 0,1,-1   / mode 0 is saving the start distribution (.start)
                    / mode 1 is saving the distribution at the current timestep
                    / mode -1 is saving the distribution at the last timestep
    If the filename is speciefied, it is used and the mode is ignored
    If no filename is specified, the mode is used to generate a generic filename
    depending on the RunNumber and the Runprefix.
    """
    # --- Create array that holds the number of particles per species --- #
    npart = numpy.zeros(top.ns,'i')
    for i in range(top.ns):
        npart[i] = getn(js=i)

    # --- Calculate beam current taking into account the scraped particles --- #
    if ibeam_old is not None and np_old is not None:
        ibeam_s = ibeam_old*npart/top.sp_fract/np_old
    else:
        ibeam_s = numpy.ones(top.ns)*largepos

    # --- Find maximum number of particles of all species --- #
    npmax = npart.max()

    # --- Set up data dict to hold particle information. --- #
    # Note: x,y,z,vx,vy,vz arrays will be of size ns, npmax. For species with
    # less particles than npmax, the empty spaces are filled with largepos (i.e.
    # 1e+36). The new current per species is calculated with the initial current
    # times the ratio of injected particles and particles crossing zpend.
    DATA = dict(M = top.aion_s,
                mode = "par",
                Q = top.zion_s,
                I = ibeam_s,
                ns = top.ns,
                np = npart,
                npmax = npmax,
                fract = top.sp_fract,
                x = numpy.ones([top.ns,npmax],'d')*largepos,
                y = numpy.ones([top.ns,npmax],'d')*largepos,
                z = numpy.ones([top.ns,npmax],'d')*top.zbeam*1000, # (mm)
                vx = numpy.ones([top.ns,npmax],'d')*largepos,
                vy = numpy.ones([top.ns,npmax],'d')*largepos,
                vz = numpy.ones([top.ns,npmax],'d')*largepos )

    for i in range(top.ns):
        ar_size = getx(js=i).size
        DATA["x"][i,:ar_size]   = getx(js=i)*1000               # (mm)
        DATA["y"][i,:ar_size]   = gety(js=i)*1000               # (mm)
        DATA["vx"][i,:ar_size]  = getvx(js=i)                   # (m/s)
        DATA["vy"][i,:ar_size]  = getvy(js=i)                   # (m/s)
        DATA["vz"][i,:ar_size]  = getvz(js=i)                   # (m/s)

    # --- Try and see if a global verbosity value has been assigned --- #
    try:
        global verbosity
        verbosity
    except: verbosity = 0

    if filename is None:
        if mode == 0:
            if verbosity == 1: print "saving initial particle distribution"
            filename = "./"+RunPrefix+RunNumber+"/"+RunPrefix+RunNumber+".start"
        elif mode == 1:
            if verbosity == 1: print "saving intermediate particle distribution"
            filename =      "./"+RunPrefix+RunNumber+"/"+RunPrefix+RunNumber+"_"+\
                                    "%1.2e"%(top.zbeam)+"m.dist"
        elif mode == -1:
            if verbosity == 1: print "saving final particle distribution"
            filename = "./"+RunPrefix+RunNumber+"/"+RunPrefix+RunNumber+".end"

    # --- Open file (write binary for better cross-platform performance) --- #
    FILE = open(filename, 'wb')
    cPickle.dump(DATA, FILE)
    FILE.close()

    return 0

def save_particles_3D(INPUT, IONS, filename = None, zcp = None):
    """
    save_particles_3D is a function to conveniently save particles at the zcp
    position.

    Note: A ZCrossingParticles object must be given to this function!

    If the filename is speciefied, it is used and the mode is ignored
    If no filename is specified, the mode is used to generate a generic filename
    depending on the RunNumber and the RunPrefix.
    """
    if filename is None:
        try: filename = "./Results/"+INPUT["runid"]+'_'+INPUT["mode"]+".end"
        except:
            print "Error in save_particles_3D: Could not determine filename!"
            return 1
    if zcp is None:
        print "Error in save_particles_3D: No ZCrossingParticles object!"
        return 1

    # --- Create array that holds the number of particles per species --- #
    npart = numpy.zeros(top.ns,'i')

    for i in range(top.ns):
        npart[i] = zcp.getn(js=i)

    # --- Calculate beam current taking into account the scraped particles --- #
    # CAVE: top.np_s has to be set manually in the main script!
    ibeam_s = IONS["current"]*npart/top.sp_fract/IONS["np"]

    # --- Find maximum number of particles of all species --- #
    npmax = npart.max()

##      fma()
##      zcp.pxxp(js=-1)
##      refresh()
##      breakpoint()
##      zcp.pyyp(js=-1)
##      refresh()
##      breakpoint()

    # --- Set up data dict to hold particle information. --- #
    # Note: x,y,z,vx,vy,vz arrays will be of size ns, npmax. For species with
    # less particles than npmax, the empty spaces are filled with largepos (i.e.
    # 1e+36). The new current per species is calculated with the initial current
    # times the ratio of injected particles and particles crossing zpend.
    # Cave: it is important that top.np_s is set manually in the main script!
    DATA = dict(    M               = top.aion_s,
                                    mode    = "par",
                                    Q               = top.zion_s,
                                    I               = ibeam_s,
                                    ns              = top.ns,
                                    np              = npart,
                                    npmax   = npmax,
                                    fract   = top.sp_fract,
                                    x               = numpy.ones([top.ns,npmax],'d')*largepos,
                                    y               = numpy.ones([top.ns,npmax],'d')*largepos,
                                    z               = numpy.ones([top.ns,npmax],'d')*zcp.zz*1000, # (mm)
                                    vx              = numpy.ones([top.ns,npmax],'d')*largepos,
                                    vy              = numpy.ones([top.ns,npmax],'d')*largepos,
                                    vz              = numpy.ones([top.ns,npmax],'d')*largepos)

    for i in range(top.ns):
        ar_size = zcp.getx(js=i).size
        DATA["x"][i,:ar_size]   = zcp.getx(js=i)*1000           # (mm)
        DATA["y"][i,:ar_size]   = zcp.gety(js=i)*1000           # (mm)
        DATA["vx"][i,:ar_size]  = zcp.getvx(js=i)                       # (m/s)
        DATA["vy"][i,:ar_size]  = zcp.getvy(js=i)               # (m/s)
        DATA["vz"][i,:ar_size]  = zcp.getvz(js=i)                       # (m/s)

    # --- Open file (write binary for better cross-platform performance) --- #
    FILE = open(filename, 'wb')
    cPickle.dump(DATA, FILE)
    FILE.close()

    return 0

def save_summary(filename = None, summary = None, results = None):
    """
    simple function that saves a text into a file
    """
    if filename is None:
        try: filename = "./"+RunPrefix+RunNumber+"/"+RunPrefix+RunNumber+".info"
        except:
            print "Error in function save_summary: Filename was None and could \
not generate generic filename. No output file written."
            return 1

    FILE = open(filename, 'w')
    FILE.write("Simulation: "+RunPrefix+RunNumber+"\n\n")
    if summary is not None:
        FILE.write(summary)
    if results is not None:
        FILE.write("\n")
        FILE.write(results)

    FILE.close()

    return 0

def save_summary_3D(INPUT, IONS, filename = None):
    """
    Write textfile with all the information of the run
    """
    if filename is None:
        filename        = './Results/'+INPUT["runid"]+'_'+INPUT["mode"]+'.inf'

    f       = open(filename,'w')
    f.write("INPUT array:\n")
    for key in INPUT.keys():
        f.write(key.ljust(20)+" = %s\n"%(str(INPUT[key])))
    f.write("\nIONS array:\n")
    f.write("number of species = %s\n"%(str(IONS["ns"])))
    f.write("\naion (amu)    zion (e)      current (A)   vbeam (m/s)\n")
    for i in range(top.ns):
        f.write(("%i"%(IONS["aion"][i])).ljust(14)+("%i"%(IONS["zion"][i])).ljust(14)+("%.4e"%(IONS["current"][i])).ljust(14)+"%.4e\n"%(IONS["vbeam"][i]))
    f.close()

    return 0

def XY_save_envelope_old(f1 = None, f2 = None):
    """
    Function to save envelope data. In order not to slow down the simulation,
    the file is not opened and closed all the time but passed from outside the
    function...
    """
    if f1 is not None and f2 is not None:
        for i in range(top.ns):         # iterate over all species
            # --- get data for saving --- #
            x               = getx(js=i)
            y               = gety(js=i)
            Np              = len(x)
            # --- Calculate the beam envelopes --- #
            # calculate mean values (mm)
            XMean   = sum(x)/Np
            YMean   = sum(y)/Np
            # (mm^2)
            XSqMean = sum(x**2)/Np
            YSqMean = sum(y**2)/Np
            # calculate standard deviations (mm)
            XSd             = sqrt(XSqMean-XMean**2)
            YSd             = sqrt(YSqMean-YMean**2)
            # calculating 1-Rms beam edges (mm)
            XMin    = XMean-2.0*XSd
            XMax    = XMean+2.0*XSd
            YMin    = YMean-2.0*YSd
            YMax    = YMean+2.0*YSd
            # --- Save emittance and envelope every timestep --- #
            if len(x) != 0:
                f1.write("%e %e %e %e %e "\
                                %(top.zbeam, XMin, XMax, YMin, YMax))
                f2.write("%e %e %e %e %e "\
                                %(top.zbeam, min(x), max(x), min(y), max(y)))
            else:
                f1.write("%e %e %e %e %e "\
                                %(top.zbeam, largepos, largepos, largepos, largepos))
                f2.write("%e %e %e %e %e "\
                                %(top.zbeam, largepos, largepos, largepos, largepos))
        f1.write("\n")
        f2.write("\n")
        return 0
    else:
        print "Error: XY_save_envelope: Not enough filenames given to function!"
        return 1

def XY_init_envelope(ns, nsteps):
    """
    Function to initialize the envelope struct and allocate necessary memory
    """
    nsteps+=1
    Envelopes = {   "mode"  : "env",                    # mode ("par" or "env")
                                    "ns"    : ns,                                   # Number of species
                                    "Nz"    : numpy.ones(ns)*nsteps,        # Number of z positions
                                    "M"             : top.aion_s,               # mass of species (amu)
                                    "Q"             : top.zion_s,               # charge state (e)
                                    "z"             : numpy.zeros([ns,nsteps]),     # z-positions (m)
                                    "xmin"  : numpy.zeros([ns,nsteps]),     # (m)
                                    "xmax"  : numpy.zeros([ns,nsteps]),     # (m)
                                    "xmean" : numpy.zeros([ns,nsteps]),     # Centroid position (m)
                                    "xrms"  : numpy.zeros([ns,nsteps]),     # 1-RMS std (m)
                                    "ymin"  : numpy.zeros([ns,nsteps]),     # (m)
                                    "ymax"  : numpy.zeros([ns,nsteps]),     # (m)
                                    "ymean" : numpy.zeros([ns,nsteps]),     # Centroid position (m)
                                    "yrms"  : numpy.zeros([ns,nsteps]),     # 1-RMS std (m)
                            }

    return Envelopes

def XY_write_envelope(Envelopes, nstep):
    """
    Function to collect envelope data during runtime
    """
    for i in range(top.ns):         # iterate over all species
        # --- get data for saving --- #
        x               = getx(js=i)
        y               = gety(js=i)
        Np              = getn(js=i)
        if Np > 0:
            # --- Calculate the beam envelopes --- #
            # calculate mean values (mm)
            XMean   = sum(x)/Np
            YMean   = sum(y)/Np
            # (mm^2)
            XSqMean = sum(x**2)/Np
            YSqMean = sum(y**2)/Np
            # calculate standard deviations (mm)
            XSd             = sqrt(XSqMean-XMean**2)
            YSd             = sqrt(YSqMean-YMean**2)
            # --- Save emittance and envelope every timestep --- #
            Envelopes["z"][i,nstep]         = top.zbeam
            Envelopes["xmin"][i,nstep]      = x.min()
            Envelopes["xmax"][i,nstep]      = x.max()
            Envelopes["xmean"][i,nstep]     = XMean
            Envelopes["xrms"][i,nstep]      = XSd
            Envelopes["ymin"][i,nstep]      = y.min()
            Envelopes["ymax"][i,nstep]      = y.max()
            Envelopes["ymean"][i,nstep]     = YMean
            Envelopes["yrms"][i,nstep]      = YSd
    return 0

def XY_save_envelope(DATA, filename = None):
    """
    Function to save the collected envelope data in a cPickle file
    """
    if filename is not None:
        FILE = open(filename, 'wb')
        cPickle.dump(DATA, FILE)
        FILE.close()
    else:
        print "Error in XY_save_envelope: no filename specified"
    return 0

def save_envelope_3D(zstart = None, zend = None, filename = None):
    """
    Function to collect emittance and envelope data for WARP RZ and 3D
    simulations
    zstart  = z position where the envelope starts
    zend    = z position where the envelope ends
    """
    if zstart is None:
        try: zstart = w3d.zmmin
        except:
            print "Error in save_envelope_3D: Could not determine zstart"
            return 1
    if zend is None:
        try: zend = w3d.zmmax
        except:
            print "Error in save_envelope_3D: Could not determine zend"
            return 1
    if filename is None:
        try: filename = "./Results/"+arraytostr(top.runid)+setup.pnumb+".env"
        except:
            print "Error in save_envelope_3D: Could not determine filename"
            return 1

    # --- set up array with z positions for gathering envelope data --- #
    zpositions      = arange(zstart, zend, w3d.dz)
    ns                      = top.ns
    nsteps          = len(zpositions)

    Envelopes = {   "mode"  : "env",                                                # mode ("par" or "env")
                                    "ns"    : ns,                                           # Number of species
                                    "Nz"    : numpy.ones(ns)*nsteps,                # Number of z positions
                                    "M"             : top.aion_s,                                   # mass of species (amu)
                                    "Q"             : top.zion_s,                                   # charge state (e)
                                    "z"             : numpy.zeros([ns,nsteps]),             # z-positions (m)
                                    "xmin"  : numpy.zeros([ns,nsteps]),             # (m)
                                    "xmax"  : numpy.zeros([ns,nsteps]),             # (m)
                                    "xmean" : numpy.zeros([ns,nsteps]),             # Centroid position (m)
                                    "xrms"  : numpy.zeros([ns,nsteps]),             # 1-RMS std (m)
                                    "ymin"  : numpy.zeros([ns,nsteps]),             # (m)
                                    "ymax"  : numpy.zeros([ns,nsteps]),             # (m)
                                    "ymean" : numpy.zeros([ns,nsteps]),             # Centroid position (m)
                                    "yrms"  : numpy.zeros([ns,nsteps]),             # 1-RMS std (m)
                            }

    # --- iterate over all z positions --- #
    for j in zpositions:
        center  = j + w3d.dz*0.5
        # --- iterate over all species --- #
        for i in range(top.ns):
            # --- get data for saving --- #
            # Note: particle data will be gathered by WARP in the window
            # [zc-wz*w3d.dz;zu+wz*w3d.dz]
            x               = getx( js=i, zc=center, wz=.5)
            y               = gety( js=i, zc=center, wz=.5)
            vx              = getvx(js=i, zc=center, wz=.5)
            vy              = getvy(js=i, zc=center, wz=.5)
            vz              = getvz(js=i, zc=center, wz=.5)
            npart   = len(x)

            # --- Calculate the emittances and beam envelopes  --- #
            if npart > 0:
                emi             = Emittances(top.aion_s[i],x,y,vx,vy,vz)
                Envelopes["z"][i,j]             = center
                Envelopes["xmin"][i,j]  = x.min()
                Envelopes["xmax"][i,j]  = x.max()
                Envelopes["xmean"][i,j] = emi["XMean"]
                Envelopes["xrms"][i,j]  = emi["XSd"]
                Envelopes["ymin"][i,j]  = y.min()
                Envelopes["ymax"][i,j]  = y.max()
                Envelopes["ymean"][i,j] = emi["YMean"]
                Envelopes["yrms"][i,j]  = emi["YSd"]

    # --- Save envelope as cPickle file --- #
    FILE = open(filename, 'wb')
    cPickle.dump(Envelopes, FILE)
    FILE.close()

    return 0

def write_set_file(filename = None):
    """
    Write a .set file containing the filenames of envelope, information
    and end distribution for easier loading into the postprocessor
    """
    if filename is not None:
        f       = open(filename, 'wb')
        f.write(os.path.splitext(os.path.split(filename)[1])[0]+".env\n")
        f.write(os.path.splitext(os.path.split(filename)[1])[0]+".inf\n")
        f.write(os.path.splitext(os.path.split(filename)[1])[0]+".end\n")
        f.close()
    return 0

# -----------------------------------------------------------------------------#
# ------------------------- POTENTIAL MAP MANIPULATION ------------------------#
# -----------------------------------------------------------------------------#
def read_VENUS_B_file(  fn = "F1mm_185_153_153_465.bin",
                                                N = numpy.array([151,151,541])):
    """
    Define a function to read in the magnetic fieldmap for the VENUS source --- #
    Note: This is preliminary, later it should be able to read in any fieldmap, but
    I still have to figure out how to set up the GUI input for that DW'11
    """
    # --- open the binary file --- #
    f = open(fn,'rb')

    # Note: First six values are x0,y0,z0, dx,dy,dz
    r0 = numpy.array(unpack('ddd', f.read(24)))
    dr = numpy.array(unpack('ddd', f.read(24)))

    # --- Allocate storage for Magnetic field data --- #
    Bx = numpy.zeros(N,'d')
    By = numpy.zeros(N,'d')
    Bz = numpy.zeros(N,'d')

    # --- Read in Dipole data from file --- #
    for i in range(N[0]):
        for j in range(N[1]):
            for k in range(N[2]):
                # --- Read in field values, change sign and convert to Tesla --- #
                Bx[i,j,k],By[i,j,k],Bz[i,j,k] = -numpy.array(unpack('ddd', f.read(24)))*1.0e-4
    f. close()

    # --- Consolidate the information in a dict --- #
    B = dict(       r0      = r0,
                            dr      = dr,
                            N       = N,
                            Bx      = Bx,
                            By      = By,
                            Bz      = Bz)
    return B

def generateELfield3D(  ELFieldV, ELFieldZStart, ELFieldZEnd, ELFieldNx,
                                                ELFieldNy, ELFieldNz, ELFieldDx, ELFieldDy,
                                                ELFieldDz, autorun=True):
    """
    Calculates a 3D meshed field from a RZ symmetric field by
    linear interpolation between 2 adjacent meshpoints. This
    can be upgraded later with a more accurate method if necessary.

    NOTE: This method is not working correctly! No idea why...DW'11

    Cave: Du to the 4-fold symmetry of the field only a quater is generated.
    """
    # --- Check input for obvious errors --- #
    if not ELFieldV.any():

        FieldData = None
        print "Error in generateELfield3D: 0-field passed to function!"
        return FieldData

    if      ELFieldV.shape[0] != ELFieldNx + 3 or\
            ELFieldV.shape[2] != ELFieldNz + 3:

        FieldData = None
        print "Error in generateELfield3D: ELFieldV shape mismatch!"
        return FieldData

    if ELFieldZStart + ELFieldDz*ELFieldNz != ELFieldZEnd:

        FieldData = None
        print "Error in generateELfield3D: zstart + nz * dz doesn't match zend!"
        return FieldData

    # Cave: Potential in RZ mode has two safety cells in "y" direction (indices
    # 0 and 2). The correct potential has the y index 1! There are also safety
    # cells in x and z: I believe it's 1 cell at low z and two cells at high z,
    # same for x.

    V2D = ELFieldV[:,1,:]

    # --- Generate matrix for holding the 3D field --- #
    # Note: Due to the RZ symmetry of the original field and WARPs capability
    # of using 4-fold symmetry in 3D simulations, only a quarter in x-y of the
    # field will be created. WARPs safety cells are omitted

    Nx = ELFieldNx
    Ny = ELFieldNx
    Nz = ELFieldNz

    V3D = numpy.zeros([Nx+3,Ny+3,Nz+3],'d')

    ELFieldXMin = 0
    ELFieldYMin = 0

    for k in range(Nz+2):
        for j in range(Ny+2):
            for i in range(Nx+2):
                r = sqrt(i*i+j*j)
                if r <= ELFieldNx+1:
                    r_min = floor(r)
                    V3D[i,j,k]      =       V2D[r_min,k] + (V2D[r_min,k] -\
                                                    V2D[r_min+1,k])*(r - r_min)


    if autorun == False:
        print "ELFieldV.shape", ELFieldV.shape
        print "V2D.shape", V2D.shape
        print "V3D.shape", V3D.shape
        breakpoint()

    FieldData = {   "V3D"   : V3D,
                                    "Nx"    : Nx,
                                    "Ny"    : Ny,
                                    "Nz"    : Nz,
                                    "Dx"    : ELFieldDx,
                                    "Dy"    : ELFieldDx,
                                    "Dz"    : ELFieldDz,
                                    "XStart": ELFieldXStart,
                                    "YStart": ELFieldYStart,
                                    "ZStart": ELFieldZStart,
                                    "ZEnd"  : ELFieldZEnd   }

    return FieldData

def join_fields(DynField3D, DynZStart, ELField3D, ELZStart, PositiveHV, ZStart):
    """
    Function takes the Einzellens and Dynamitron fields, performs
    some checks and joins them into one large 3D field map for WARP.

    PositiveHV is the HV potential of the terminal of the accelerator.

    Cave: This function requires that the potential maps are both 3D in 4-fold
              symmetry and the output field will also be only a quarter of the full
              potential map.
              Furthermore the cell sizes must match!
    """

    # Perform a check on the cell sizes:
    # Cave: Due to the fact that Nx,Ny,Nz have to be even numbers,
    # Dx,Dy,Dz might not be *exactly* the same even if they were
    # set the same in the scripts because they are recalculated.
    # Hence the rounding...
    if      round(DynField3D["Dx"], 5) != round(ELField3D["Dx"], 5) or\
    round(DynField3D["Dy"], 5) != round(ELField3D["Dy"], 5) or\
    round(DynField3D["Dz"], 5) != round(ELField3D["Dz"], 5):
        print   "Error in join_fields: Cell sizes of the fields don't match!"
        print   "ELField  (dx, dy, dz): ",\
                        ELField3D["Dx"],ELField3D["Dx"],ELField3D["Dx"]
        print   "DynField (dx, dy, dz): ",\
                        DynField3D["Dx"], DynField3D["Dx"], DynField3D["Dx"]
        print   "Aborting Simulation"
        return None

    # Define some variables for convenience
    Dx = DynField3D["Dx"]
    Dy = DynField3D["Dy"]
    Dz = DynField3D["Dz"]

    # Calculate size of the map
    Nx = max([DynField3D["Nx"],ELField3D["Nx"]])
    Ny = max([DynField3D["Ny"],ELField3D["Ny"]])
    NDynStart       = int(floor((DynZStart-ELZStart)/Dz))
    Nz = NDynStart + DynField3D["Nz"] - 1
    # Note: The pot maps have three safety cells in z direction (one at 0, 2 at
    #               max z)

    # Create empty map on HV potential
    # Cave: PositiveHV is given in kV and has to be made float by putting 1000.0
    JoinedField = numpy.ones([Nx+3,Ny+3,Nz+3],'d')*PositiveHV*1000.0

    if autorun == False:
        print "ELField3D[Nx]", ELField3D["Nx"]
        print "ELField3D[Ny]", ELField3D["Ny"]
        print "ELField3D[Nz]", ELField3D["Nz"]
        print ELField3D["V3D"].shape

        print "DynField3D[Nx]", DynField3D["Nx"]
        print "DynField3D[Ny]", DynField3D["Ny"]
        print "DynField3D[Nz]", DynField3D["Nz"]
        print DynField3D["V3D"].shape

        print "JoinedField.shape", JoinedField.shape
        print "NDynStart/End", NDynStart, NDynStart+DynField3D["Nz"]
        breakpoint()

    # Add fields to map
    JoinedField[0:ELField3D["Nx"]+3, 0:ELField3D["Ny"]+3, 0:ELField3D["Nz"]+1]=\
    ELField3D["V3D"][:,:,0:-2]+PositiveHV*1000.0
    # Note: One has to be careful not to add the safety cells too!
    JoinedField[0:DynField3D["Nx"]+3,0:DynField3D["Ny"]+3, NDynStart:NDynStart+\
    DynField3D["Nz"]+2] = DynField3D["V3D"][:,:,1:]

    Full3DPot = {   "V3D"   : JoinedField,
                "Nx"        : Nx,
                                    "Ny"    : Ny,
                                    "Nz"    : Nz,
                                    "Dx"    : Dx,
                                    "Dy"    : Dy,
                                    "Dz"    : Dz,
                                    "XStart": ELFieldXStart,
                                    "YStart": ELFieldYStart }

    return Full3DPot

# -----------------------------------------------------------------------------#
# ------------------------------ INPUT FUNCTIONS ------------------------------#
# -----------------------------------------------------------------------------#

def boltzmann(np,boltzmannE):
    """
    generate random velocity with Boltzmann distribution
    """

    Velo = dict(vx=zeros(np, 'd'),
                vy=zeros(np, 'd'),
                vz=zeros(np, 'd'),
                eV=boltzmannE)

    for i in range(np):

        phi=random()*2.0*pi
        theta=(random())*pi/2.0001      # pi/2.0 so z always positive
        extraE=2.0*echarge*boltzmannE/top.aion_s[w3d.inj_js]/amu
        temp=sqrt(-extraE*log(random()))
        Velo["vy"][i]=temp*sin(phi)
        Velo["vx"][i]=temp*cos(phi)
        Velo["vz"][i]=sqrt(-extraE*log(random())*cos(theta))
##              print abs(Velo["vx"]/Velo["vz"]).mean()
##              print abs(Velo["vy"]/Velo["vz"]).mean()
##              breakpoint()
    return Velo

def load_particles(filename):
    """
    Function to load particles from a .end or .dist file
    """
    # --- Load the particle distribution from a .dist or .end file --- #
    if os.path.splitext(filename)[1] not in [".dist", ".end", ".start", ".mid"]:
        print "Error in load_particles: Wrong extension, must be .end or .dist"
        raise SystemExit
    try:
        FILE = open(filename, 'rb')
    except IOError:
        print "Couldn't open the particle distribution file. Aborting simulation"
        raise SystemExit
    try:
        DISTRIBUTION = cPickle.load(FILE)
    except EOFError:
        print "Couldn't unpickle the particle distribution file, maybe not a cPickle file? Aborting simulation"
        raise SystemExit
    FILE.close()

    return DISTRIBUTION

def load_SuSI_solenoid(Solenoid_Path, autorun = False):
    """
    Loads the solenoid data on axis into a WARP array
    """
    with open(Solenoid_Path, 'rb') as FILE:

        B_DATA = []

        for line in FILE.readlines():

            values = line.split(',')

            # Note: position in file is in mm, field in T
            B_DATA.append((float(values[0])*1.0e-3, float(values[1])))

    B_DATA = np.array(B_DATA, dtype = [('z', float), ('Bz', float)])

    nz = len(B_DATA['z'])
    dz = B_DATA['z'][1] - B_DATA['z'][0]
    Bmax_ind = np.where(abs(B_DATA['Bz']) == max(abs(B_DATA['Bz'])))[0][0]

    if not autorun:
        print "Setting source solenoid field"
        print "Field starts at %.2f mm"%(B_DATA['z'][0]*1000.0)
        print "dz = %.2f mm, nz = %i, Bz max = %.2f T at z = %.2f mm"%(
            dz*1000.0, nz, B_DATA['Bz'][Bmax_ind], B_DATA['z'][Bmax_ind]*1000.0)

        breakpoint()

    top.mmltid[0] = 1
    top.mmltzs[0] = 0.0 # Start of field is at simulation zero reference
    top.mmltze[0] = (nz-1)*dz # end of field (m)
    top.mmltap[0] = .15
    top.mmltsf[0] = 0.0
    top.mmltsc[0] = 1.0

    top.nmmltsets = 1
    top.nmsmult = 1
    top.nzmmltmax = nz
    gchange("Mult_data")
    top.nzmmlt[0] = nz
    top.dzmmlt[0] = dz

    top.msmmlt[0:nz,0,0] = B_DATA['Bz']

    return 0

def createRZbeam(NumPart, ns, radius):
    """
    Creates an RZ beam with similar rms values as the triangular initial
    distribution in ECR extraction simulations.
    """

    particles = {"np": np.ones(ns, 'int') * NumPart,
                 "x": np.zeros([ns, NumPart],'d'),
                 "y": np.zeros([ns, NumPart],'d'),
                 "vx": np.zeros([ns, NumPart],'d'),
                 "vy": np.zeros([ns, NumPart],'d'),
                 "vz": np.zeros([ns, NumPart],'d')}

    for i in range(ns):

        r = radius * np.sqrt(np.random.random(NumPart))
        theta = 2.0 * np.pi * np.random.random(NumPart)

        particles["x"][i][:] = r * np.cos(theta)
        particles["y"][i][:] = r * np.sin(theta)
        particles["vz"][i][:] = np.ones(NumPart, 'd') * IONS["vbeam"][i]

    return particles

def rzInjection():
    """
    This function will be called once per species during injection to create the
    particles to be injected. The current species number is stored in w3d.inj_js
    The particles are created on the surface of the
    source (hence there are no z values specified). The code gets the
    data from the xt, yt etc arrays. Note that these arrays may be zeroed
    out afterward since they are used as temporary arrays for other code chunks
    """

    # --- Set number of particles in species --- #
    NumPart = particles["np"][w3d.inj_js]
    w3d.npgrp = NumPart

    # --- (Re-)allocate space for particles --- #
    gchange('Setpwork3d')

    # --- Load particles into WARP's temporary arrays --- #
    w3d.xt[:] = particles["x"][w3d.inj_js][:NumPart]
    w3d.yt[:] = particles["y"][w3d.inj_js][:NumPart]
    w3d.uxt[:] = particles["vx"][w3d.inj_js][:NumPart]
    w3d.uyt[:] = particles["vx"][w3d.inj_js][:NumPart]
    w3d.uzt[:] = particles["vx"][w3d.inj_js][:NumPart]
    print
    print "Currents:", top.ibeam_s, "A"
    print "Total current:", sum(top.ibeam_s), "A"
    print "sum of top.curr:", sum(top.curr)
    return 0

def myinjection():
    """
    This function will be called once per species during injection to create the particles
    to be injected. The current species number is stored in w3d.inj_js
    The particles are created on the surface of the
    source (hence there are no z values specified). The code gets the
    data from the xt, yt etc arrays. Note that these arrays may be zeroed
    out afterward since they are used as temporary arrays for other code chunks
    """
    mode = INPUT["mode"]

    if mode == "triangle":

        np_temp = DISTRIBUTION["np"][w3d.inj_js]

        for i in range(DISTRIBUTION["np"][w3d.inj_js]):

            if sqrt(DISTRIBUTION["x"][w3d.inj_js,i]**2+DISTRIBUTION["y"][w3d.inj_js,i]**2) >= INPUT["sourceradius"]*1000:

                np_temp = np_temp - 1

        w3d.npgrp = np_temp

    else: w3d.npgrp = int(NumPart)

    gchange("Setpwork3d")

    xt  = numpy.zeros(shape(w3d.xt),'d')
    yt  = numpy.zeros(shape(w3d.xt),'d')
    uxt = numpy.zeros(shape(w3d.uxt),'d')
    uyt = numpy.zeros(shape(w3d.uxt),'d')

    # --- rotation of the triangle due to misalignment of sextupole magnet --- #
    # --- turn off for SuSI
##    rotang=(-90-16.0)*pi/180.0

    if mode == "triangle":

        print 'reading in particle data for species '+str(w3d.inj_js+1)

        index = 0

        for i in range(DISTRIBUTION["np"][w3d.inj_js]):

            if sqrt(DISTRIBUTION["x"][w3d.inj_js, i]**2+DISTRIBUTION["y"][w3d.inj_js, i]**2) < INPUT["sourceradius"]*1000:

                xt[index]  = DISTRIBUTION["x"][w3d.inj_js,i]*1.0e-3
                yt[index]  = DISTRIBUTION["y"][w3d.inj_js,i]*1.0e-3
                uxt[index] = DISTRIBUTION["vx"][w3d.inj_js,i]
                uyt[index] = DISTRIBUTION["vy"][w3d.inj_js,i]

                index += 1

##              #Cave: Temporary for testing (uncomment uxt and uyt above)
##              velo = boltzmann(np_temp,INPUT["iontemp"])
##              for i in range(np_temp):
##                      uxt[i]  = velo["vx"][i]
##                      uyt[i]  = velo["vy"][i]
    elif mode in ["circle", "IGUN"]:

        # --- Round distributions for testing --- #
        velo = boltzmann(NumPart,INPUT["iontemp"])

        for i in range(NumPart):

            radius  = sqrt(random())*INPUT["sourceradius"]
            angle   = random()*2*pi
            xt[i]   = radius*numpy.cos(angle)
            yt[i]   = radius*numpy.sin(angle)
            uxt[i]  = velo["vx"][i]
            uyt[i]  = velo["vy"][i]

    # --- Set the WARP arrays for the particles --- #
    # Note: rotation doesn't change round beam...
    w3d.xt = xt; w3d.yt = yt
    w3d.uxt = uxt; w3d.uyt = uyt
##    w3d.xt  =        xt*numpy.cos(rotang)   -  yt*numpy.sin(rotang)
##    w3d.yt  =        xt*numpy.sin(rotang)   +  yt*numpy.cos(rotang)
##    w3d.uxt =       uxt*numpy.sin(rotang)   + uyt*numpy.cos(rotang)
##    w3d.uyt =       uxt*numpy.sin(rotang)   + uyt*numpy.cos(rotang)
    # --- Set all ions of one species to the same zvelocity (ion acoustic speed)
    # Note: top.vzinject will be added to uzt when w3d.l_inj_user_particles_v       = False,
    # so it needs to be set appropriately. Note further that if top.vzinject is zero (and
    # inject is set to constant current) it will be set to top.vbeam during the generate,
    # and so will need to be reset.
    if mode == "IGUN":

        w3d.uzt = 0

    else:

        w3d.uzt = top.vbeam_s[w3d.inj_js]

    return 0

def ReadIGUN(filename=None, nparticles=5000):
    """
    Function to read in an IGUN .TRJ file and calculate particle distribution
    data from it.
    """
    if filename is not None:
    # extract basename (incl. path) and extension
        path,ext=os.path.splitext(filename)
        if ext in [".trj",".TRJ",".Trj"]:
            try: f = open(filename, 'r')
            except IOError:
                print "Could not open file!"
                return -1
        else:
            print "Wrong file extension! Must be '.TRJ'!"
            return -1
    else: return -1

    run_label = f.readline()
    data = f.readlines()
    f.close()
    Nlines = len(data)-1

    # Allocate some arrays:
    ray             = 0
    group           = numpy.zeros(Nlines, 'd')
    charge          = numpy.zeros(Nlines, 'd')
    mass            = numpy.zeros(Nlines, 'd')
    r                       = numpy.zeros(Nlines, 'd')
    z_temp          = 0
    rp                      = numpy.zeros(Nlines, 'd')
    raycurrent      = numpy.zeros(Nlines, 'd')
    energy          = numpy.zeros(Nlines, 'd')
    phip            = numpy.zeros(Nlines, 'd')
    phi         = numpy.zeros(Nlines, 'd')

    for i in range(Nlines):
        ray, group[i], charge[i], mass[i], r[i], z_temp, energy[i], rp[i],\
        raycurrent[i], phip[i], phi[i]  = data[i].split()

    ns = int(max(group))     # number of species
    pps = int(Nlines/ns)     # particles per species

    # Reshape the arrays according to the number of species (flat --> matrices):
    group           = numpy.reshape(group, (ns,pps))
    charge          = numpy.reshape(charge, (ns,pps))
    mass            = numpy.reshape(mass, (ns,pps))
    r                       = numpy.reshape(r, (ns,pps))
    rp                      = numpy.reshape(rp, (ns,pps))
    raycurrent      = abs(numpy.reshape(raycurrent, (ns,pps)))
    energy          = numpy.reshape(energy, (ns,pps))
    phip            = numpy.reshape(phip, (ns,pps))
    phi         = numpy.reshape(phi, (ns,pps))

    x                       = numpy.zeros((ns,nparticles), 'd')
    y                       = numpy.zeros((ns,nparticles), 'd')
    xp                      = numpy.zeros((ns,nparticles), 'd')
    yp                      = numpy.zeros((ns,nparticles), 'd')
    z                       = numpy.zeros((ns,nparticles), 'd')
    vx                      = numpy.zeros((ns,nparticles), 'd')
    vy                      = numpy.zeros((ns,nparticles), 'd')
    vz          = numpy.zeros((ns,nparticles), 'd')
    currentsum  = numpy.zeros(ns,'d')

    for j in range(ns):         # process each species
        currentsum[j] = sum(raycurrent[j,:])
        cumulative = zeros(pps+1,'d')
        for i in range(pps):
            cumulative[i+1]=cumulative[i]+raycurrent[j,i]/currentsum[j]

        for i in range(nparticles):
            probability = numpy.random.random_sample()  # get random number
            jmin=0
            jmid=int(pps/2)
            jmax=pps

            for dummy in range(40):
                if cumulative[jmin] <= probability and probability <=\
                                                                                                cumulative[jmid]:
                    if jmin+1 == jmid:
                        jmid=jmin
                        break
                    jmax = jmid
                    jmid = int((jmin+jmax)/2)
                elif cumulative[jmid] <= probability and probability <=\
                                                                                                cumulative[jmax]:
                    if jmid+1 == jmax: break
                    jmin = jmid
                    jmid = int((jmin+jmax)/2)
                else:
                    print "%s: probability %g of out boundaries cumulative[%d]\
                               = %g - cumulative[%d] = %g\n"\
                               %(os.path.split(sys.argv[0])[1],probability,jmin,\
                               cumulative[jmin],jmax,cumulative[jmax])

            jmid-=1

            theta = 2.0*numpy.pi*numpy.random.random_sample()
            velocity =      sqrt(2.0*energy[j,jmid]*1.602176487e-19*charge[j,0]/\
                                    (mass[j,0]*1.66053886e-27))

            x[j,i]  = r[j,jmid]*numpy.cos(theta)                                                    # (mm)
            y[j,i]  = r[j,jmid]*numpy.sin(theta)                                                    # (mm)
            z[j,i]  = z_temp                                                                                        # (mm)
            xp[j,i] = (rp[j,jmid]*numpy.cos(theta)-phip[j,jmid]*\
                                                                                                    numpy.sin(theta))       # (rad)
            yp[j,i] = (rp[j,jmid]*numpy.sin(theta)+phip[j,jmid]*\
                                                                                                    numpy.cos(theta))       # (rad)
            vz[j,i] = velocity/sqrt(xp[j,i]**2+yp[j,i]**2+1)                        # (m/s)
            vx[j,i] = xp[j,i]*vz[j,i]                                                                       # (m/s)
            vy[j,i] = yp[j,i]*vz[j,i]                                                                       # (m/s)

    # Calculate some handy additional output values
    vzmean  = vz.mean(axis=1)       # Calculate mean vz for each species (m/s)
    # Note: The mean velocity is multiplied by the sqrt of the charge state,
    # for some reason IGUN does not take the charge state of multiple species
    # into account ??? DW'11
    # CAVE: What happens if the charge state in INPUT5 is not 1?

    xmax    = x.max(axis=1)
    ymax    = y.max(axis=1)

    xenv = numpy.zeros(ns, 'd')
    yenv = numpy.zeros(ns, 'd')
    for i in range(ns):
        maxid = where(x[i,:]==xmax[i])
        xenv[i] = xp[i,maxid[0]]        # rad
        maxid = where(y[i,:]==ymax[i])
        yenv[i] = yp[i,maxid[0]]        # rad

    results = {     "ns"                    : ns,
                                    "pps"                   : pps,
                                    "M"                             : mass[:,0],
                                    "Q"                             : charge[:,0],
                                    "totalCurrent"  : abs(currentsum*6.2832),
                                    "x"                             : x,
                                    "y"                             : y,
                                    "z"                             : z,
                                    "xp"                    : xp*1000, # mrad
                                    "yp"                    : yp*1000, # mrad
                                    "vx"                    : vx,
                                    "vy"                    : vy,
                                    "vz"                    : vz,
                                    "vzmean"                : vzmean,
                                    "xmax"                  : xmax,
                                    "ymax"                  : ymax,
                                    "xenv"                  : xenv,
                                    "yenv"                  : yenv   }

    return results

# -----------------------------------------------------------------------------#
# ------------------------------- MISCELLANEOUS -------------------------------#
# -----------------------------------------------------------------------------#
def plotfields(no_of_steps=1):
    """
    Function for printing out the conductors and fields + particles
    """
    # make graph every "no_of_steps" (defaults to 10) timesteps
    if top.it%no_of_steps==0:       # top.it is the global timestep counter
        fma()
        conductors.draw()
        pfzx(plotsg=0,cond=0, phicolor=blue, titles=0)
        #if 1: limits(0, 0.15, -0.025, 0.025)    # good overview over the extraction simulation
        if 1: limits(0, 0.15, -0.032, 0.032)    # ANP: See all of electrode        
        if 0: limits(-0.009506,0.05,-0.0075,0.0075) # magnify the plasma region
        gunppzx(js=-1, titles=0)
        #ppzx(js=-1,titles=0)                        # plot particles
        pfzx(plotsg=0,cond=0, rhocolor=red, phicolor=green, contours = linspace(24000, w3d.plasmapotential, 10), titles=0)  # plot meniscus
        refresh()

    return 0

def plotfields_old(no_of_steps=1):
    """
    Function for printing out the conductors and fields + particles
    """
    # make graph every "no_of_steps" (defaults to 10) timesteps
    if top.it%no_of_steps==0:       # top.it is the global timestep counter
        fma()
        conductors.draw(titles=0)
        pfzx(plotsg=0,cond=0, phicolor=blue, titles=0)
        if 1: limits(0, 0.15, -0.025, 0.025)    # good overview over the extraction simulation
        if 0: limits(-0.009506,0.05,-0.0075,0.0075) # magnify the plasma region
        gunppzx(js=-1, titles=0)
        #ppzx(js=-1,titles=0)                        # plot particles
        pfzx(plotsg=0,cond=0, rhocolor=red, phicolor=green, contours = linspace(24000, w3d.plasmapotential, 10), titles=0)  # plot meniscus
        refresh()

    return 0

def plot_xy():
    """
    Function that plots the xy particle distribution after the first timestep
    and waits for the user
    """
    print "SW:", top.pgroup.sw
    if top.pgroup.zp.mean() > -0.0081:
        fma()
        ppxxp(js=-1)
        refresh()
        if not autorun: breakpoint()
        ppyyp(js=-1)
        refresh()
        if not autorun: breakpoint()
        uninstallafterstep(plot_xy)

    return 0

def goodnn(nn):
    """
    Function that makes the gridcell number divisible by 2 for FFT (Fast Fourier
    Transformation
    """
    if nn%2 == 0: return nn
    else:         return nn + 1

def setplasmarho():
    """
    This function fills in the charge density in the rest of the plasma
    volume. It copies the density from near the edge of the beam and
    replicates it out to the edge of the grid.
    """
    iz = int((zsrcmin[0]-w3d.zmmin)/w3d.dz)
    ix = int(top.a0/w3d.dx) - 1
    frz.basegrid.rho[ix:,:iz+1] = frz.basegrid.rho[ix,:iz+1]

    return 0

def redirect_stdout():
    """
    Redirect stdout to stringIO object to read in number of v-cycles
    """
    sys.stdout = mystdout

    return 0

def get_vcycles():
    """
    Function to catch the stdout line from the fieldsolver, search for a
    significant word ('v-cycles') and output the error and number of iterations
    for further use...
    """
    sys.stdout = old_stdout
    global mystdout
    output = mystdout.getvalue()
    mystdout.seek(0)
    mystdout.truncate(0)
    print output
    global mgiters
    global mgerror
    v_index = output.find("v-cycles")
    mgiters = int(output[v_index-5:v_index-1])
    mgerror = float(output[v_index-20:v_index-9])

    return 0

def breakpoint(msg = "continue? Enter/Esc..."):
    """
    Convenience function to pause the code and ask for continuation or abortion
    Cave: Does not work on Windows Systems!
    """
    print msg

    fd = sys.stdin.fileno()

    oldterm = termios.tcgetattr(fd)
    newattr = termios.tcgetattr(fd)
    newattr[3] = newattr[3] & ~termios.ICANON & ~termios.ECHO
    termios.tcsetattr(fd, termios.TCSANOW, newattr)

    oldflags = fcntl.fcntl(fd, fcntl.F_GETFL)
    fcntl.fcntl(fd, fcntl.F_SETFL, oldflags | os.O_NONBLOCK)

    try:
        while 1:
            try:
                c = sys.stdin.read(1)
                if c == '\x1b': quit()
                elif c== '\n': break
            except IOError: pass
    finally:
        termios.tcsetattr(fd, termios.TCSAFLUSH, oldterm)
        fcntl.fcntl(fd, fcntl.F_SETFL, oldflags)

def EX_make_output_folders(CWD):
    """
    Makes the necessary output folders for extraction simulations
    """
    RESULTS_PATH = os.path.join(CWD,"Results")
    if not os.path.isdir(RESULTS_PATH):
        os.mkdir(RESULTS_PATH)

    RESULTS_PATH = os.path.join(RESULTS_PATH, "Extraction")
    if not os.path.isdir(RESULTS_PATH):
        os.mkdir(RESULTS_PATH)

    RELAXED_PATH    = os.path.join(CWD,"RelaxedPotentials")
    if not os.path.isdir(RELAXED_PATH):
        os.mkdir(RELAXED_PATH)

    return CWD, RESULTS_PATH, RELAXED_PATH

def SuSI_EX_Output_Folders(CWD):
    """
    Makes the necessary output folders for extraction simulations
    """
    RESULTS_PATH = os.path.join(CWD,"Results")
    if not os.path.isdir(RESULTS_PATH):
        os.mkdir(RESULTS_PATH)

    RELAXED_PATH = os.path.join(CWD,"RelaxedPotentials")
    if not os.path.isdir(RELAXED_PATH):
        os.mkdir(RELAXED_PATH)

    return CWD, RESULTS_PATH, RELAXED_PATH

def KillParticles():
    """
    Very simple particle scraper by Dr. Albe Lemut 2010

    Usage:  Use 'installparticlescraper(KillParticles) in the WARP XY Slice
                    code. During the step by step simulation, one can change the value
                    of SimulationTubeRadius.
                    If no SimulationTubeRadius has been specified, the Scraper will use
                    w3d.xmmax as outermost boundary to scrape the particles.
    """
    try:
        global SimulationTubeRadius
        SimulationTubeRadius
    except:
        SimulationTubeRadius = w3d.xmmax

    r = sqrt(top.pgroup.xp**2+top.pgroup.yp**2)
    top.pgroup.gaminv[r>SimulationTubeRadius] = 0

    return 1

def KillParticles2():
    """
    Very simple particle scraper by Dr. Albe Lemut 2010
    Version 2 by D Winklehner 2013

    Usage:  Use 'installparticlescraper(KillParticles2) in the WARP XY Slice
            code. During the step by step simulation, one can change the
            value stored in Limits = [xmin, ymin, xmax, ymax, rmax].
            If no Limits have been specified, the Scraper will use
            w3d.xmmax as rmax to scrape the particles.
    """
    # Limits consists of the 5 values [xmin, xmax, ymin, ymax and r],
    # if either of them is 'None', we assume there is only a round pipe
    # If Limits is not specified, we use w3d.xmmax as pipe radius
    try:

        global Limits
        Limits

    except:

        Limits = [None, None, None, None, w3d.xmmax]

    if None in Limits:

        r = sqrt(top.pgroup.xp**2+top.pgroup.yp**2)
        top.pgroup.gaminv[r>Limits[4]] = 0

    else:

        xmin, xmax, ymin, ymax, rmax = Limits

        top.pgroup.gaminv[top.pgroup.xp>xmax] = 0
        top.pgroup.gaminv[top.pgroup.xp<xmin] = 0
        top.pgroup.gaminv[top.pgroup.yp>ymax] = 0
        top.pgroup.gaminv[top.pgroup.yp<ymin] = 0

    return 0

def Kicker(xangle, yangle):
    """
    """
    x_angle = xangle/180.0*np.pi
    y_angle = yangle/180.0*np.pi

    # Apply horizontal kick:
    uzp_new = top.pgroup.uzp*np.cos(x_angle) - \
              top.pgroup.uxp*np.sin(x_angle)

    top.pgroup.uxp = top.pgroup.uzp*np.sin(x_angle) + \
                     top.pgroup.uxp*np.cos(x_angle)

    top.pgroup.uzp = uzp_new

    #Apply vertical kick
    uzp_new = top.pgroup.uzp*np.cos(y_angle) - \
              top.pgroup.uyp*np.sin(y_angle)

    top.pgroup.uyp = top.pgroup.uzp*np.sin(y_angle) + \
                     top.pgroup.uyp*np.cos(y_angle)

    top.pgroup.uzp = uzp_new

    return 0

def Emittance(M,x,y,vx,vy,vz,printscreen=0):
    """
    Albe's original Emittance calculation function
    """
    x               = x*1000.0 # (mm)
    y               = y*1000.0 # (mm)

    # calculate x and y prime (mrad)
    xp              = vx/vz*1000.0
    yp              = vy/vz*1000.0

    # get number of alive particles
    Np              = len(x)

    # calculate mean velocity module (m/s)
    V               = sqrt(vx**2+vy**2+vz**2)
    VMean   = sum(V)/Np
    VSd             = sqrt(sum((V-VMean)**2)/Np)

    # light speed (m/s)
    c               = 2.99792458e8

    # calculate beta
    beta    = VMean/c

    # calculate gamma
    gamma   = 1.0/sqrt(1.0-beta*beta)

    # calculate beta*gamma for normalized emittance
    BetaGamma       = beta*gamma

    # amu to kg
    amu             = 1.66053886e-27

    # calculate mean kinetic energy (keV)
    T               = 0.5*M*amu*V**2/1.6022e-16
    TMean   = sum(T)/Np
    TSd             = sqrt(sum((T-TMean)**2)/Np)

    # calculate mean values

    # (mm)
    XMean=sum(x)/Np
    YMean=sum(y)/Np

    # (mrad)
    XPMean=sum(xp)/Np
    YPMean=sum(yp)/Np

    # (mm*mrad)
    XXPMean=sum(x*xp)/Np
    YYPMean=sum(y*yp)/Np

    # calculate mean square values

    # (mm^2)
    XSqMean=sum(x**2)/Np
    YSqMean=sum(y**2)/Np

    # (mrad^2)
    XPSqMean=sum(xp**2)/Np
    YPSqMean=sum(yp**2)/Np

    # (mm^2*mrad^2)
    XXPSqMean=sum((x*xp)**2)/Np
    YYPSqMean=sum((y*yp)**2)/Np

    # calculate standard deviations

    # (mm)
    XSd=sqrt(XSqMean-XMean**2)
    YSd=sqrt(YSqMean-YMean**2)

    # (mrad)
    XPSd=sqrt(XPSqMean-XPMean**2)
    YPSd=sqrt(YPSqMean-YPMean**2)

    # (mm*mrad)
    AXXPSd=sqrt(XXPSqMean-XXPMean**2)
    AYYPSd=sqrt(YYPSqMean-YYPMean**2)

    # (mm*mrad)
    XXPSd=sum((x-XMean)*(xp-XPMean))/Np
    YYPSd=sum((y-YMean)*(yp-YPMean))/Np

    # calculating 2-Rms beam edges (mm)
    XMin=XMean-2.0*XSd
    XMax=XMean+2.0*XSd

    YMin=YMean-2.0*YSd
    YMax=YMean+2.0*YSd

    XPMin=XPMean-2.0*XPSd
    XPMax=XPMean+2.0*XPSd

    YPMin=YPMean-2.0*YPSd
    YPMax=YPMean+2.0*YPSd

    # calculating 2-Rms beam edges prime (mm)
    XMinP=XPMean-2.0*XXPSd/XSd
    XMaxP=XPMean+2.0*XXPSd/XSd

    YMinP=YPMean-2.0*YYPSd/YSd
    YMaxP=YPMean+2.0*YYPSd/YSd

    # calculating 1-Rms Emittances (pi*mm*mrad)
    XERms=sqrt((XSd*XPSd)**2-XXPSd**2)/pi
    YERms=sqrt((YSd*YPSd)**2-YYPSd**2)/pi

    # calculating 1-Rms Normalized Emittances (pi*mm*mrad)
    XENRms=XERms*BetaGamma
    YENRms=YERms*BetaGamma

    # calculating 1-Rms Normalized Emittances (pi*mm*mrad)
    XENRms=XERms*BetaGamma
    YENRms=YERms*BetaGamma

    # calculating 4-Rms Emittances (pi*mm*mrad)
    XE4Rms=4.0*XERms
    YE4Rms=4.0*YERms

    # calculating 4-Rms Normalized Emittances (pi*mm*mrad)
    XEN4Rms=XE4Rms*BetaGamma
    YEN4Rms=YE4Rms*BetaGamma

    # Calculating 4Rms Courant Snyder parameters

    XBeta=XMax**2/XE4Rms
    XGamma=XMaxP**2/XE4Rms
    XAlphaSign=XXPMean/abs(XXPMean)
    if XBeta*XGamma < 1.0:
        XAlpha=0.0
    else:
        XAlpha=-XAlphaSign*sqrt(XBeta*XGamma-1.0)

    YBeta=YMax**2/XE4Rms
    YGamma=YMaxP**2/XE4Rms
    YAlphaSign=YYPMean/abs(YYPMean)
    if YBeta*YGamma < 1.0:
        YAlpha=0.0
    else:
        YAlpha=-YAlphaSign*sqrt(YBeta*YGamma-1.0)

    xl=x-XMean
    yl=y-YMean
    rl=sqrt(xl**2+yl**2)
    cosl=xl/rl
    sinl=yl/rl
    re=sqrt((2*XSd*cosl)**2+(2*YSd*sinl)**2)
    NpEnv=0
    for i in range(Np):
        if rl[i] <= re[i]:
            NpEnv=NpEnv+1

    RMax=max(sqrt(x**2+y**2))
    RPMax=max(sqrt(xp**2+yp**2))

    # printing results on screen if required
    if printscreen == 1:

        print "    Np      =",Np
        print "    VMean   =",VMean,"(m/s) VSd =",VSd,"(m/s) Bg =",BetaGamma
        print "    TMean   =",TMean,"(keV) TSd =",TSd,"(keV)"
        print "    XMean   =",XMean,"(mm)", "XSd =",XSd,"(mm)"
        print "    YMean   =",YMean,"(mm)", "YSd =",YSd,"(mm)"
        print "    XPMean  =",XPMean,"(mrad)", "XPSd =",XPSd,"(mrad)"
        print "    YPMean  =",YPMean,"(mrad)", "YPSd =",YPSd,"(mrad)"
        print "    XXPMean =",XXPMean,"(mm*mrad)", "XXPSd =",XXPSd,"(mm*mrad)"
        print "    YYPMean =",YYPMean,"(mm*mrad)", "YYPSd =",YYPSd,"(mm*mrad)"
        #print "    XXPMean =",XXPMean,"(mm*mrad)", "XXPSd =",XXPSd,"(mm*mrad) AXXPSd =",AXXPSd,"(mm*mrad)"
        #print "    YYPMean =",YYPMean,"(mm*mrad)", "YYPSd =",YYPSd,"(mm*mrad) AYYPSd =",AYYPSd,"(mm*mrad)"
        print "    XMin    =",XMin,"(mm)", "XMax =",XMax,"(mm)"
        print "    YMin    =",YMin,"(mm)", "YMax =",YMax,"(mm)"
        print "    XPMin   =",XPMin,"(mm)", "XPMax =",XPMax,"(mm)"
        print "    YPMin   =",YPMin,"(mm)", "YPMax =",YPMax,"(mm)"
        print "    XMinP   =",XMinP,"(mrad) XMaxP =",XMaxP,"(mrad)"
        print "    YMinP   =",YMinP,"(mrad) YMaxP =",YMaxP,"(mrad)"
        print "    XERms   =",XERms, "(pi*mm*mrad) XENRms =", XENRms,"(pi*mm*mrad)"
        print "    YERms   =",YERms, "(pi*mm*mrad) YENRms =", YENRms,"(pi*mm*mrad)"
        print "    XE4Rms  =",XE4Rms, "(pi*mm*mrad) XEN4Rms =", XEN4Rms,"(pi*mm*mrad)"
        print "    YE4Rms  =",YE4Rms, "(pi*mm*mrad) YEN4Rms =", YEN4Rms,"(pi*mm*mrad)"
        print "    XGamma  =",XGamma,"(mrad/mm) XAlpha =", XAlpha,"(Sign =","%+.0f"%(XAlphaSign)+") XBeta =", XBeta,"(mm/mrad)"
        print "    YGamma  =",YGamma,"(mrad/mm) YAlpha =", YAlpha,"(Sign =","%+.0f"%(YAlphaSign)+") YBeta =", YBeta,"(mm/mrad)"
        print "    NpEnv   =",NpEnv
        print "    RMax    =",RMax,"(mm)"
        print "    RPMax   =",RPMax,"(mrad)"
        print

    return ( Np,TMean,XMin,XMean,XMax,YMin,YMean,YMax,XPMin,XMinP,XPMean,XMaxP,
                     XPMax,YPMin,YMinP,YPMean,YMaxP,YPMax,XERms,YERms,XENRms,YENRms,
                     NpEnv,RMax,RPMax )


def Emittances(M,x,y,vx,vy,vz, output=False):
    """
    Helper file for calculating the normalized RMS emittance of an
    ion beam from from the particles' positions and velocities.

    Function takes arrays M, x, y, z, vx, vy, vz as input and returns
    emi_rms_norm_xx', _yy' (soon also _xy', _yx') plus the beam edges for each
    species.

    Daniel Winklehner LBL 2010
    dwinklehner@lbl.gov

    Input arrays must be numpy arrays!
    """
    # m to mm
    x=x*1000.0 # (mm)
    y=y*1000.0 # (mm)

    # calculate x and y prime (mrad)
    xp=vx/vz*1000.0
    yp=vy/vz*1000.0

    # get number of alive particles
    Np=len(x)

    # calculate mean velocity module (m/s)
    V=sqrt(vx**2+vy**2+vz**2)
    VMean=sum(V)/Np
    VSd=sqrt(sum((V-VMean)**2)/Np)

    # light speed (m/s)
    c=2.99792458e8

    # calculate beta
    beta=VMean/c

    # calculate gamma
    gamma=1.0/sqrt(1.0-beta*beta)

    # calculate beta*gamma for normalized emittance
    BetaGamma=beta*gamma

    # amu to kg
    amu=1.66053886e-27

    # calculate mean kinetic energy (keV)
    T=0.5*M*amu*V**2/1.6022e-16
    TMean=sum(T)/Np
    TSd=sqrt(sum((T-TMean)**2)/Np)

    # calculate mean values

    # (mm)
    XMean=sum(x)/Np
    YMean=sum(y)/Np

    # (mrad)
    XPMean=sum(xp)/Np
    YPMean=sum(yp)/Np

    # (mm*mrad)
    XXPMean=sum(x*xp)/Np
    YYPMean=sum(y*yp)/Np

    XYPMean=sum(x*yp)/Np
    YXPMean=sum(y*xp)/Np

    # calculate mean square values

    # (mm^2)
    XSqMean=sum(x**2)/Np
    YSqMean=sum(y**2)/Np

    # (mrad^2)
    XPSqMean=sum(xp**2)/Np
    YPSqMean=sum(yp**2)/Np

    # (mm^2*mrad^2)
    XXPSqMean=sum((x*xp)**2)/Np
    YYPSqMean=sum((y*yp)**2)/Np
    XYPSqMean=sum((x*yp)**2)/Np
    YXPSqMean=sum((y*xp)**2)/Np

    # calculate standard deviations

    # (mm)
    XSd=sqrt(XSqMean-XMean**2)
    YSd=sqrt(YSqMean-YMean**2)

    # (mrad)
    XPSd=sqrt(XPSqMean-XPMean**2)
    YPSd=sqrt(YPSqMean-YPMean**2)

    # (mm*mrad)
    AXXPSd=sqrt(XXPSqMean-XXPMean**2)
    AYYPSd=sqrt(YYPSqMean-YYPMean**2)
    AXYPSd=sqrt(XYPSqMean-XYPMean**2)
    AYXPSd=sqrt(YXPSqMean-YXPMean**2)

    # (mm*mrad)
    XXPSd=sum((x-XMean)*(xp-XPMean))/Np
    YYPSd=sum((y-YMean)*(yp-YPMean))/Np
    XYPSd=sum((x-XMean)*(yp-YPMean))/Np
    YXPSd=sum((y-YMean)*(xp-XPMean))/Np

    # Calculating 1-rms beam edges (mm)
    XMin1rms=XMean-XSd
    XMax1rms=XMean+XSd

    YMin1rms=YMean-YSd
    YMax1rms=YMean+YSd

    # calculating 2-Rms beam edges (mm)
    XMin2rms=XMean-2.0*XSd
    XMax2rms=XMean+2.0*XSd

    YMin2rms=YMean-2.0*YSd
    YMax2rms=YMean+2.0*YSd

    # Calculating max envelope (mm)
    XMinFull=x.min()
    XMaxFull=x.max()

    YMinFull=y.min()
    YMaxFull=y.max()

    # calculating 2-Rms beam edges prime (mm)
    XMinP=XPMean-2.0*XXPSd/XSd
    XMaxP=XPMean+2.0*XXPSd/XSd

    YMinP=YPMean-2.0*YYPSd/YSd
    YMaxP=YPMean+2.0*YYPSd/YSd

    # calculating 1-Rms Emittances (mm*mrad)
    # --- WARNING: There has been some controversy between Albe and myself about
    # --- dividing by pi here or not! DW'10
    XERms=sqrt((XSd*XPSd)**2-XXPSd**2)
    YERms=sqrt((YSd*YPSd)**2-YYPSd**2)
    XYPRms=sqrt((XSd*YPSd)**2-XYPSd**2)
    YXPRms=sqrt((YSd*XPSd)**2-YXPSd**2)

    # calculating 1-Rms Normalized Emittances (mm*mrad)
    XENRms=XERms*BetaGamma
    YENRms=YERms*BetaGamma
    XYPENRms=XYPRms*BetaGamma
    YXPENRms=YXPRms*BetaGamma

    # calculating 4-Rms Emittances (pi*mm*mrad)
    XE4Rms=4.0*XERms
    YE4Rms=4.0*YERms

    # calculating 4-Rms Normalized Emittances (pi*mm*mrad)
    XEN4Rms=XE4Rms*BetaGamma
    YEN4Rms=YE4Rms*BetaGamma

    # Calculate Twiss Params (DW) ThinkAboutThis! DW'11
    Xbeta=(XSd**2)/XERms
    Xgamma=(XPSd**2)/XERms
    if XXPSd < 0:
        Xalpha = sqrt(Xbeta*Xgamma-1)
    else:
        Xalpha = -sqrt(Xbeta*Xgamma-1)

    Ybeta=(YSd**2)/YERms
    Ygamma=(YPSd**2)/YERms
    if YYPSd < 0:
        Yalpha = sqrt(Ybeta*Ygamma-1)
    else:
        Yalpha = -sqrt(Ybeta*Ygamma-1)

    # Calculating 4-Rms Courant Snyder parameters
##      XBeta=XMax**2/XE4Rms
##      XGamma=XMaxP**2/XE4Rms
##      XAlphaSign=XXPMean/abs(XXPMean)
##      if XBeta*XGamma < 1.0:
##              XAlpha=0.0
##      else:
##              XAlpha=-XAlphaSign*sqrt(XBeta*XGamma-1.0)
##
##      YBeta=YMax**2/XE4Rms
##      YGamma=YMaxP**2/XE4Rms
##      YAlphaSign=YYPMean/abs(YYPMean)
##      if YBeta*YGamma < 1.0:
##          YAlpha=0.0
##      else:
##              YAlpha=-YAlphaSign*sqrt(YBeta*YGamma-1.0)

    # Calculating the max emittance by using the rms Twiss parameters and
    # the maximum values of x, x', y and y'

    EXXPMax = Xgamma*max(x)*max(x)+2*Xalpha*max(x)*max(xp)+Xbeta*max(xp)*max(xp)
    EYYPMax = Ygamma*max(y)*max(y)+2*Yalpha*max(y)*max(yp)+Ybeta*max(yp)*max(yp)

    Summary="""
--- RMS Emittances, Twiss Parameters and Beamsizes ---
(Emittances are 1-RMS)

Np             = %i
v_mean         = %1.5e m/s
E_kin_mean     = %1.5e keV

e_rms_xxp      = %5.4f mm-mrad
e_rms_yyp      = %5.4f mm-mrad
e_rms_xyp      = %5.4f mm-mrad
e_rms_yxp      = %5.4f mm-mrad

e_rms_norm_xxp = %5.4f mm-mrad
e_rms_norm_yyp = %5.4f mm-mrad
e_rms_norm_xyp = %5.4f mm-mrad
e_rms_norm_yxp = %5.4f mm-mrad

e_full_xxp     = %5.4f pi-mm-mrad
e_full_yyp     = %5.4f pi-mm-mrad
(calculated from ellipse that includes all particles and the Twiss parameters)

alpha_x        = %5.4f
beta_x         = %5.4f mm/mrad
gamma_x        = %5.4f mrad/mm

alpha_y        = %5.4f
beta_y         = %5.4f mm/mrad
gamma_y        = %5.4f mrad/mm

dia_x_1rms     = %5.4f mm
dia_y_1rms     = %5.4f mm

dia_x_2rms     = %5.4f mm
dia_y_2rms     = %5.4f mm

dia_x_full     = %5.4f mm
dia_y_full     = %5.4f mm
    """%(   Np, VMean, TMean,
                    XERms, YERms, XYPRms, YXPRms,
                    XENRms, YENRms, XYPENRms, YXPENRms,
                    EXXPMax, EYYPMax,
                    Xalpha, Xbeta, Xgamma,
                    Yalpha, Ybeta, Ygamma,
                    XMax1rms-XMin1rms, YMax1rms-YMin1rms,
                    XMax2rms-XMin2rms, YMax2rms-YMin2rms,
                    XMaxFull-XMinFull, YMaxFull-YMinFull    )

    if output == True:
        print Summary

    Results = dict( Np = Np,
                    VMean           = VMean,
                    TMean           = TMean,
                    EXXPRMS         = XERms,
                    EYYPRMS         = YERms,
                    EXYPRMS         = XYPRms,
                    EYXPRMS         = YXPRms,
                    EXXPNRMS        = XENRms,
                    EYYPNRMS        = YENRms,
                    EXYPNRMS        = XYPENRms,
                    EYXPNRMS        = YXPENRms,
                    XMin1rms        = XMin1rms,
                    XMax1rms        = XMax1rms,
                    YMin1rms        = YMin1rms,
                    YMax1rms        = YMax1rms,
                    XMin2rms        = XMin2rms,
                    XMax2rms        = XMax2rms,
                    YMin2rms        = YMin2rms,
                    YMax2rms        = YMax2rms,
                    XMinFull        = XMinFull,
                    XMaxFull        = XMaxFull,
                    YMinFull        = YMinFull,
                    YMaxFull        = YMaxFull,
                    XAlpha          = Xalpha,
                    XBeta           = Xbeta,
                    XGamma          = Xgamma,
                    YAlpha          = Yalpha,
                    YBeta           = Ybeta,
                    YGamma          = Ygamma,
                    XMean           = XMean,
                    XPMean          = XPMean,
                    YMean           = YMean,
                    YPMean          = YPMean,
                    EXXPMax = EXXPMax,
                    EYYPMax = EYYPMax,
                    Summary = Summary,
                    XSd = XSd,
                    YSd = YSd   )

##      # --- Create a structured array for better postprocessing --- #
##      # types of the entries
##      dtype   = [     ('name', str), ('value', float), ('unit', str)]
##
##      # actual entries
##      values  =       [       ('# of particles', Np, ""),
##                                      ('mean velocity', VMean, "m/s"),
##                                      ('mean energy', TMean, "keV"),
##                                      ('x-dia 1-rms', XMax1rms-XMin1rms, 'mm'),
##                                      ('y-dia 1-rms', YMax1rms-YMin1rms, 'mm'),
##                                      ('x-dia 2-rms', XMax2rms-XMin2rms, 'mm'),
##                                      ('y-dia 2-rms', YMax2rms-YMin2rms, 'mm'),
##                                      ('x-dia max', XMaxFull-XMinFull, 'mm'),
##                                      ('x-dia max', YMaxFull-YMinFull, 'mm'),
##                                      ('x centroid offset', XMean, 'mm'),
##                                      ('y centroid offset', YMean, 'mm'),
##                                      ('xp centroid offset', XPMean, 'mm'),
##                                      ('yp centroid offset', YPMean, 'mm'),
##                                      ('e xxp 1-rms', XERms, 'mm-mrad'),
##                                      ('e yyp 1-rms', YERms, 'mm-mrad'),
##                                      ('e xyp 1-rms', XYPRms, 'mm-mrad'),
##                                      ('e yxp 1-rms', YXPRms, 'mm-mrad'),
##                                      ('e norm xxp 1-rms', XENRms, 'mm-mrad'),
##                                      ('e norm yyp 1-rms', YENRms, 'mm-mrad'),
##                                      ('e norm xyp 1-rms', XYPENRms, 'mm-mrad'),
##                                      ('e norm yxp 1-rms', YXPENRms, 'mm-mrad'),
##                                      ('4e xxp 1-rms', 4.0*XERms, 'mm-mrad'),
##                                      ('4e yyp 1-rms', 4.0*YERms, 'mm-mrad'),
##                                      ('4e xyp 1-rms', 4.0*XYPERms, 'mm-mrad'),
##                                      ('4e yxp 1-rms', 4.0*YXPERms, 'mm-mrad'),
##                                      ('4e norm xxp 1-rms', 4.0*XENRms, 'mm-mrad'),
##                                      ('4e norm yyp 1-rms', 4.0*YENRms, 'mm-mrad'),
##                                      ('4e norm xyp 1-rms', 4.0*XYPENRms, 'mm-mrad'),
##                                      ('4e norm yxp 1-rms', 4.0*YXPENRms, 'mm-mrad'),
##                              ]
##
##      # put together in structured array
##      emi_array = numpy.array(values, dtype = dtype)

    return Results

def coulomblog(ne, q, vb):
    """
    returns the coulomb logarithm
    """
    echarge  = 1.602e-19       # unit charge in C
    amu      = 1.6605e-27      # atomic mass unit in kg
    epsilon0 = 8.854187817e-12 # vacuum permittivity in F/m
    clight   = 3e8             # speed of light in vacuum in m/s
    me       = 9.10938291e-31  # electron mass in kg
    kB       = 1.3806488e-23   # Boltzmann constant m^2*kg/s^2/K
    Tr       = 300             # Room temperature in K (~300 K)
    hbar     = 1.05457173e-34  # PLanck constant

    b_90 = q*echarge**2/4/np.pi/epsilon0/me/vb**2

    om_pe = np.sqrt(ne*echarge**2/(me*epsilon0))

    b_max = vb/om_pe

    L = 4*np.pi*np.log(b_max/b_90)

    return L

def neutralization_multi(n0, I, sigma_e, sigma_i, q, m, vb, r0):
    """
    Calculate the neutralization according to Gabovich's formula
    extended for multiple species
    this assumes that all inputs are numpy arrays of the same length
    (r0 will be averaged)

    """
    echarge  = 1.602e-19       # unit charge in C
    amu      = 1.6605e-27      # atomic mass unit in kg
    epsilon0 = 8.854187817e-12 # vacuum permittivity in F/m
    clight   = 3e8             # speed of light in vacuum in m/s
    me       = 9.10938291e-31  # electron mass in kg
    kB       = 1.3806488e-23   # Boltzmann constant m^2*kg/s^2/K
    Tr       = 300             # Room temperature in K (~300 K)
    hbar     = 1.05457173e-34  # PLanck constant
    PHIi     = 15.51
    Ti = 10 # secondary ion temperature in eV
    mi = 28  # ~Mass of N2 in amu
    vi = np.sqrt(echarge*Ti*2/(mi*amu))

    V = 24000

    r0 = np.average(r0)

##    print "I", I

    delta_phi_full = np.sum(I/(4*np.pi*vb*epsilon0))

##    print "delta_phi_full", delta_phi_full

    nb = I/(echarge*q*vb*np.pi*r0**2)
    ni = r0*n0*nb*vb*sigma_i/(2*vi)
    ne = np.sum(q*nb + ni)

##    print "ne", ne

    L = coulomblog(ne, q, vb)

##    print "L", L

    chi = np.sqrt(ne/n0*np.sum(I*m*amu*L)/np.sum(I*sigma_e/q)*echarge**2/(4*np.pi*epsilon0)**2*PHIi*3/me/V)/delta_phi_full

    neut = (1 + 0.5*chi**2 - 0.5*chi*np.sqrt(4+chi**2))

    return neut

class FourierFit:
    """
    CAVE: NOT WORKING PROPERLY AT THE MOMENT!
    """

    def __init__(self, order = 0):
        self.order  =   order
        self.a      =   numpy.ones(order+1,'d')
        self.b      =   numpy.ones(order+1,'d')
        self.w      =   1.0
        self.limits =   numpy.zeros(2,'d')

    def set_order(self, order = 0):
        self.order  =   order
        self.a      =       numpy.ones(order+1,'d')
        self.b      =   numpy.ones(order+1,'d')
        self.w      =   1.0

        return 0

    def f(self, x):

        if x.min() < self.limits[0] or x.max() > self.limits[1]:

            return None

        else:

            f = numpy.zeros(x.shape,'d') + self.a[0]

            for i in range(self.order):

                f       += self.a[i+1]*numpy.cos((i+1)*x*self.w)
                f       += self.b[i+1]*numpy.sin((i+1)*x*self.w)

            return f

    def get_limits(self):
        """
        Returns the stored limits of the good fit region
        """

        return self.limits

    def set_limits(self, zlim):
        """
        zlim...tuple containing (zmin, zmax)
        """

        self.limits = numpy.array(zlim,'d')

        return 0

    def f_eval(self, x, a, b, w):
        """
        """
        f = a[0]

        for i in range(self.order):

            f       += a[i+1]*numpy.cos((i+1)*x*w)
            f       += b[i+1]*numpy.sin((i+1)*x*w)

        return f

    def initial_guess(self):
        """
        """

        params = []
        params.extend(self.a)
        params.extend(self.b)
        params.extend([self.w])

        return params

    def residuals(self, p, y, x):
        """
        p = [a,b,w]
        """

        err = y - self.f_eval(  x,
                                p[:self.order+1],
                                p[self.order+1:2*self.order+2],
                                p[-1]   )

        return err

    def least_squares_fit(self, y, x):
        """
        p = [a,b,w]
        """

        plsq    = leastsq(      self.residuals,
                                self.initial_guess(),
                                args    = (y, x),
                                maxfev  = 10000         )

        for element in plsq:
            print element

        p   = plsq[0]   # The values of the optimized parameters are stored here

        self.a  = p[:self.order+1]
        self.b  = p[self.order+1:2*self.order+2]
        self.w  = p[-1]

        return 0

    def set_DIANA_field(self):
        """
        sets the fit coefficients for the DIANA field calculated with Pandira
        """
        self.order = 8
        self.a = numpy.array([291,2568,2289,748.9,-22.28,-119.8,-51.04,-12.28,-3.352],'d')
        self.b = numpy.array([0,-1192,402.2,662.7,226.3,-9.588,-12.4,17.86,13.94],'d')
        self.w = 0.03085
        self.limits = numpy.array([-11.5, 160.0],'d')

        return 0

    def set_LEDA_field(self):
        """
        sets the fit coefficients for the DIANA field calculated with Pandira
        """
        self.order = 8
        self.a = numpy.array([698.5,813.1,169.4,12.8,1.246,10.21,6.552, 2.176,0.2586],'d')
        self.b = numpy.array([0,-446.9,-261.2,-80.43,-7.952,1.641,-0.7795,-1.556,-0.6301],'d')
        self.w = 0.008503
        self.limits = numpy.array([-300.0, 300.0],'d')

        return 0

