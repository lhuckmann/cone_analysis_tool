import numpy as np
import cone_analysis_tool 
import time

cat=cone_analysis_tool.cone_analysis_tool()

global_start=time.time()
print('START CONE ANALYSIS\n')
print('atom ID', 'cone angle')

#---- set up cell ---
cell_const=14.6488
xyz = cat.read_xyz('geom_opt.xyz')
cell = [[cell_const, 0, 0], [0.0, cell_const, 0], [0.0, 0.0, cell_const]]

r_cut = 2.25 # cut-off for the neighbor list


#---- construct NL ---
start=time.time()
ids, nl = cat.nl(xyz, xyz, r_cut, cell)
nl_time=time.time()-start

#---- cone analysis ---
TIME=[]
cone=[]
for ii in ids:
    start=time.time()
    _, _, _, angle = cat.best_spot(xyz[ii], np.array(xyz[nl[ii]]), cell, 1, 2000)
    print(ii, np.round(angle,2))
    cone.append(angle)
    TIME.append(time.time()-start)

cone=np.array(cone)

#---- print results ---
print(40*'+')
print('NL time:', np.round(nl_time,5), '\nsampling time:', np.round(np.sum(TIME),5), '\nsampling time/atom:', np.round(np.mean(TIME),5), '\n')
print('Total runtime:', np.round(time.time()-global_start,5))

np.savetxt('cone_angles.dat', np.array([np.arange(0,len(cone),1), cone]).T, header='idx\t theta [deg]', fmt=['%u', '%1.6f'])
