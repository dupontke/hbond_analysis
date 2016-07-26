#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# USAGE:

# PREAMBLE:

from plotting_functions import *
from sel_list import *

dat = sys.argv[1]  
system = sys.argv[2]

nSel = len(sel)

# ----------------------------------------
# MAIN PROGRAM:

datalist = np.loadtxt(dat)

nSteps = len(datalist[:])               # for 2D use the following: nSteps = len(datalist[:,0]             
print 'Number of selections: %d, number of steps: %d' %(nSel,nSteps)

time = np.zeros(nSteps)

for i in range(nSteps):
	time[i] = i*0.002		# units of time in ns; each frame is separated by 0.002 ns 

for i in range(nSel):
	selection = sel[i][0]

	scat_hist(time[:],datalist[:],'k','Time (ns)','RMSD','%02d.%s' %(i,selection),'%s' %(system),yunits='$\AA$')


# for the last for loop datalist[:,i]
