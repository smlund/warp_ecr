#!/usr/bin/env python

import os
import numpy as npy
import matplotlib.pyplot as plt
import cPickle

def load_particles(filename):
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

initParticlePath='./SourceSuSI_v1_000_10.end'

distri=load_particles(initParticlePath)

# Change value of index to see other charge states
index=0

x=distri['x'][0,0:distri['np'][0]]
y=distri['y'][0,0:distri['np'][0]]

#plt.plot(x,y,'ro',markersize=0.4)
plt.scatter(x,y,s=0.5,color='red',edgecolor='')
plt.title('Ions of '+'{:.3f}'.format(distri['M'][index])+' AMU and '+'{:d}'.format(int(distri['Q'][index]))+'+ charge state')
plt.xlabel('x (mm)')
plt.ylabel('y (mm)')
plt.xlim(-60,60)
plt.ylim(-60,60)
plt.gca().set_aspect('equal',adjustable='box')
plt.show()