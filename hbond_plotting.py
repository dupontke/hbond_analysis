#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# USAGE:

# ./hbond_plotting.py data_file system_descriptor

# PREAMBLE:

from plotting_functions import *
from sel_list import *

dat = sys.argv[1]  
system = sys.argv[2]

nSel = len(sel)

# ----------------------------------------
# MAIN PROGRAM:

# Load in .dat_file
datalist = np.loadtxt(dat)
hbonds = len(datalist[:])
print 'Number of selections: %d, Total hydrogen bonds found: %d' %(nSel,hbonds)


# Run the plotting function for Hydrogen Bond Analysis.
# Descriptions for scat_hist(xdata, ydata, color, x_axis, y_axis, system, analysis, **kwargs)
# NOTE: [:,0] refers to [row,column]. [:] will loop through all the rows. [0] loop through the first column in all the rows. 
for i in range(nSel):
	selection = sel[i][0]
	scat_hist(datalist[:,0],datalist[:,1],'k','Time (ns)','Hydrogen Bond Distance','%02d.%s' %(i,selection),'%s' %(system),yunits='$\AA$',x_lim=(0,40))
