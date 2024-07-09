import numpy as np
import cone_analysis_tool 
import time

cat=cone_analysis_tool.cone_analysis_tool()

global_start=time.time()
print('START CONE ANALYSIS\n')
print('atom ID', 'cone angle')

#---- set up cell ---
cell_const=14.6488
elem, xyz = cat.read_xyz('geom_opt.xyz')
cell = [[cell_const, 0, 0], [0.0, cell_const, 0], [0.0, 0.0, cell_const]]

r_cut = 2.25 # cut-off for the neighbor list

dict_cut={'Si': {'Si': 2.7, 'O': 2.00, 'N': 2.25},
          'N': {'Si': 2.25, 'O': 2.25, 'N': 2.25},
          'O': {'Si': 2.0, 'O': 2.25, 'N': 2.25}}

#---- construct NL ---
start=time.time()
idxs, nl = cat.nl(xyz, xyz, elem, elem, dict_cut, cell)
nl_time=time.time()-start

#---- cone analysis ---
TIME=[]
cone=[]
for ii in idxs:
    start=time.time()
    _, _, _, angle = cat.best_spot(xyz[ii], np.array(xyz[nl[ii]]), cell, 1, 2000)
    print(ii, np.round(angle,2))
    cone.append(angle)
    TIME.append(time.time()-start)

cone=np.array(cone)
cn=cat.get_CN(nl)
print(cn)
#---- print results ---
print(40*'+')
print('NL time:', np.round(nl_time,5), '\nsampling time:', np.round(np.sum(TIME),5), '\nsampling time/atom:', np.round(np.mean(TIME),5), '\n')
print('Total runtime:', np.round(time.time()-global_start,5))

np.savetxt('cone_angles.dat', np.array([np.arange(0,len(cone),1), cone, cn]).T, header='idx\t theta [deg]\t CN', fmt=['%u', '%1.6f', '%u'])
