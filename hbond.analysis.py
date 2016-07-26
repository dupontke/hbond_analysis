#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:

# ./hbond.analysis.py pdb_file trajectory_location start_traj end_traj system_descriptor

# ----------------------------------------
# PREAMBLE:

import numpy as np
import sys
import os
import MDAnalysis
import MDAnalysis.core.AtomGroup
import MDAnalysis.analysis.hbonds
import MDAnalysis.coordinates.base
from sel_list import *


# VARIABLE DECLARATION:

pdb = sys.argv[1]                   # point to a pdb or prmtop or psf file (untested for both prmtop and psf files)
traj_loc = sys.argv[2]              # point to the location of the trajectory files
start_traj = int(sys.argv[3])
end_traj = int(sys.argv[4])
system = sys.argv[5]

# SELECTIONS FOR HYDROGEN BOND ANALYSIS
#selection1_type = 'both'
#detect_hydrogens = 'distance'
#start = 'None'
#stop = 'None'
#step = 'None'

nSel = len(sel)

flush = sys.stdout.flush

# ----------------------------------------
# SUBROUTINES:

def ffprint(string):
    print '%s' %(string)
    flush()

def summary(nSteps):
    sum_file = open('%s.hbond.summary' %(system,sel[i][0]),'w')
    sum_file.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
    sum_file.write('To recreate this analysis, run this line in terminal:\n')
    for i in range(len(sys.argv)):
        sun_file.write('%s ' %(sys.argv[i]))
    sum_file.write('\n\n')
    sum_file.write('Progress output is written to:\n')
    sum_file.write('\nTotal number of steps analyzed: %d\n' %(nSteps))
    sum_file.write('\nAtom selections analyzed:\n')
    for i in range(nSel):
        sum_file.write('%02d  %s  %s\n' %(i,sel[i][0],sel[i][1],sel[i][2]))
    sum_file.close()


# ----------------------------------------
# MAIN:
#

# INITIALIZING UNIVERSE

u = MDAnalysis.Universe(pdb)
out_list = []
h_list = []
table_list = []
for i in range(nSel):
    h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(u, selection1 = sel[i][1], selection2 = sel[i][2], selection1_type='both', update_selection1=False, update_selection2=False, detect_hydrogens='distance', start=None, stop=None, step=None, distance=3.0, angle=120.0, donors=None, acceptors=None)
    h_list.append(h)

# BEGINNING TO ANALYZE TRAJECTORIES
nSteps = 0
while start_traj <= end_traj:
    ffprint('Loading trajectory %s' %(start_traj))
    u.load_new('%sproduction.%s/production.%s.dcd' %(traj_loc,start_traj,start_traj))
    nSteps += len(u.trajectory)
##    t0 = u.trajectory.start_timestep / u.trajectory.skip_timestep * u.trajectory.dt
##    t = [t0 + (ts.frame-1) * u.trajectory.dt for ts in u.trajectory]
    t0 = u.trajectory.start_timestep / u.trajectory.skip_timestep * 0.002
    t = [t0 + (ts.frame-1) * 0.002 for ts in u.trajectory]
#    print t[:6]
#
##    frame_numbers = [ts.frame for ts in u.trajectory]
##    print frame_numbers[:]
#    sys.exit()
    for i in range(nSel):
        out1 = open('%s.%s.results.dat' %(system,sel[i][0]), 'a')
        h_list[i].run()
        htimeseries = h_list[i].timeseries     
#        print htimeseries
#        print len(h.timeseries)
#        print h_list[i].timesteps
#        bonded_hs = h_list[i]._get_bonded_hydrogens()
#        print bonded_hs
##       htype = h_list[i].timesteps_by_type()
##        print htype
#        sys.exit()
#	print len(htimeseries)
	for j in range(len(u.trajectory)):
#        for j in range(1000):
#            print 'frame',j
##		print htimeseries[0] 		### the h-bond data for frame i... 
            if len(h_list[i].timeseries[j]) == 0:
                out1.write('%.6f    %.9f    %.9f\n' %(t[j],0.0,0.0))
#                print t[j], 0.00, 0.00
            else: 
                for k in range(len(htimeseries[j])):
##			print htimeseries[0][0]
                #for l in range(len(htimeseries[j][k])):
#                    print t[j], htimeseries[j][k][-2], htimeseries[j][k][-1]
                    out1.write('%.6f    %.9f    %.9f\n' %(t[j],htimeseries[j][k][-2],htimeseries[j][k][-1]))
        out1.close()
        
#	ffprint('Finished analyzing trajectory %02d\n' %(start_traj))
    start_traj += 1

#ffprint('Analyzed %d steps.' %(nSteps))

#summary(nSteps)

