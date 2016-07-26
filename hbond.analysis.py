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
    h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(u,selection1 = sel[i][1],selection2 = sel[i][2],selection1_type='both',detect_hydrogens='distance',start=None,stop=2,step=None, distance=3.0, angle=120.0)
    h_list.append(h)


# BEGINNING TO ANALYZE TRAJECTORIES
nSteps = 0
while start_traj <= end_traj:
    ffprint('Loading trajectory %s' %(start_traj))
    u.load_new('%sproduction.%s/production.%s.dcd' %(traj_loc,start_traj,start_traj))
    nSteps += len(u.trajectory)
    for i in range(nSel):
        h_list[i].run()
        htimeseries = h_list[i].timeseries
        
        print h_list[i].timesteps
        htype = h_list[i].timesteps_by_type()
        print htype
        sys.exit()
	print len(htimeseries)
	for j in range(len(u.trajectory)):
#		print htimeseries[0] 		### the h-bond data for frame i... 
		for k in range(len(htimeseries[i])):
#			print htimeseries[0][0]
			for l in range(len(htimeseries[j][k])):
				print htimeseries[0][0][-2], htimeseries[0][0][-1]
				sys.exit()

#        hgentable = h_list[i].generate_table()
#        table = h_list[i].table
#        table_list.append(table)
#        out1 = open('%s.%s.table.dat' %(system,sel[i][0]), 'a')
#        out1.write('%s' %(table))
#        out1.close()

        # extract number of frames and cooresponding hydrogen bond distances for the trajectories
#        out2 = open('' %(system,))
#        frames = len(htimeseries)
#        ffprint('Frames = %s\n' %(frames))
        
	ffprint('Finished analyzing trajectory %02d\n' %(start_traj))
    start_traj += 1

ffprint('Analyzed %d steps.' %(nSteps))

summary(nSteps)

