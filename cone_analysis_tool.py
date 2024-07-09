"""
_________________________________________________________

	CONE ANALYSIS TOOL
_________________________________________________________
"""
import numpy as np
import random

class cone_analysis_tool:
    def __init__(self):
        print('Welcome to the cone analysis tool')
    #--------------------------------
    #   MINIMUM IMAGE VECTOR (MIC)
    #
    # Requires:
    #   xi      | [array], shape (3,) coordinates of atom i
    #   xj      | [array], shape (3,) coordinates of atom j
    #   cell    | [array], shape (3,3) cell matrix of shape
    #           |   cell=[[ x,  0,  0],
    #                     [ xy, y,  0],
    #                     [ xz, yz, z]]
    #
    # Returns:
    #   minimum image distance, minimum image vector  
    # 
    # Notes:
    #   Cell follows the lammps semantics, see 
    #   https://docs.lammps.org/Howto_triclinic.html
    #
    
    def mic(self, x_i, x_j, cell):
        r_ij = x_j - x_i
        for ii in range(2,-1,-1):
            if r_ij[ii] > cell[ii][ii]/2:
                r_ij = r_ij-cell[ii]
            elif r_ij[ii] < -cell[ii][ii]/2:
                r_ij = r_ij+cell[ii]
            else:
                pass
        return np.sqrt(np.dot(r_ij, r_ij)), r_ij
    
    #--------------------------------
    #   PAIR DISTANCE MATRIX (PDM)
    #
    # Requires:
    #   xyz1 | array/list (3,N) coordinates of the particles of interest
    #   xyz2 | array/list (3,M) coordinates of the reference particles.
    #
    # Returns:
    #   pair distance matrix (N,M)  
    #
    # Notes:
    #   Uses two distinct coord vectors. Computes the neighborlist for the atoms
    #   in xyz1 with respect to the ones in xyz2. Allows to get the PDM i.e. for 
    #   subsets of an atomistic system.
    #   If xyz1==xyz2, the computation is accelerated by computing only the upper
    #   half of the pdm.   
    #
    def pdm(self, xyz1, xyz2, cell):
        if np.array_equal(xyz1,xyz2):
            R_ABS_MAT=np.zeros((len(xyz1),len(xyz1)))
            for ii in range(len(xyz1)):
                for jj in range(len(xyz1)):
                    if ii < jj:
                        r_ij,_ = self.mic(xyz1[ii], xyz2[jj], cell)
                        R_ABS_MAT[ii,jj]=r_ij
            R_ABS_MAT = R_ABS_MAT + R_ABS_MAT.T - np.diag(np.diag(R_ABS_MAT))
        else:
            R_ABS_MAT=np.zeros((len(xyz1),len(xyz2)))
            for ii in range(len(xyz1)):
                for jj in range(len(xyz2)):
                    r_ij,_ = self.mic(xyz1[ii], xyz2[jj], cell)
                    R_ABS_MAT[ii,jj]=r_ij
        return R_ABS_MAT
     
    #--------------------------------
    #   NEIGHBOR LIST (NL)
    #
    # Requires:
    #   xyz1      | array/list (3,N) coordinates of the particles of interest
    #   xyz2      | array/list (3,M) coordinates of the target particles.
    #   r_cut     | float () cut-off radius for 1st coordination shell
    #   cell      | [array], shape (3,3) cell matrix of shape. See MIC.
    #   R_ABS_MAT | array/list (N,M) if already given: The PDM. (optional)
    # 
    # Returns:
    #   atom indices, neighbor indices (both zero bound)
    #
    # Notes:
    #   Uses two distinct coord vectors. Computes the neighborlist for the atoms
    #   in xyz1 with respect to the ones in xyz2. Allows to get the NL i.e. for 
    #   subsets of an atomistic system.
    #
    #   Can use external PDM.
    #
    def nl(self, xyz1, xyz2, r_cut, cell, R_ABS_MAT=np.array([])):
        IDX=[]
        NL=[]
        if R_ABS_MAT.any():
            for ii in range(len(R_ABS_MAT)):
                IDX.append(ii) 
                neighs = np.where((R_ABS_MAT[ii] < r_cut) & (R_ABS_MAT[ii] > 0))[0]
                NL.append(neighs)
        else:
            R_ABS_MAT = self.pdm(xyz1, xyz2, cell)
            for ii in range(len(R_ABS_MAT)):
                IDX.append(ii) 
                neighs = np.where((R_ABS_MAT[ii] < r_cut) & (R_ABS_MAT[ii] > 0))[0]
                NL.append(neighs)
        return IDX, NL
    
    #--------------------------------
    #   READ XYZ
    #
    # Requires:
    #   filename | (string) self explaining
    # 
    # Returns:
    #   x, y, z
    # Notes:
    #   only single frame, no trajectory
    #
    def read_xyz(self,filename):
        return np.loadtxt(filename, skiprows=2, usecols=(1,2,3))
    
    
    #--------------------------------
    #   SPHERE POINT
    #
    # Requires:
    #   xyz0 | [array] (3,) center of the sphere
    #   r    | [float]  radius of the sphere
    # 
    # Returns:
    #   [x, y, z] of the sphere point
    # Notes:
    #   --
    #
    def sphere_point(self,xyz0, r):
        x=random.uniform(-r,r)
        ylim=np.sqrt(r**2-x**2)
        y=random.uniform(-ylim, ylim)
        z=np.sqrt(r**2-x**2-y**2)*random.choice((-1,1))
        return [xyz0[0]+x,xyz0[1]+y,xyz0[2]+z]
    
    #--------------------------------
    #   ANGLE
    #
    # Requires:
    #   xyz1    | [array] (3,) coordinates of the center atom
    #   xyz2    | [array] (3,) coordinates of the neighbor 1
    #   xyz3    | [array] (3,) coordinates of the neighbor 2
    #   cell    | [array] (3,3) cell matrix of shape. See MIC.
    #
    # Returns:
    #   angle (float) in degrees
    #
    # Notes:
    #   --
    #
    def angle(self,xyz1, xyz2, xyz3, cell):
        _, r12=self.mic(xyz2,xyz1,cell)
        _, r13=self.mic(xyz3,xyz1,cell)
        return np.rad2deg(np.arccos(np.dot(r12, r13)/(np.sqrt(np.dot(r12, r12))*np.sqrt(np.dot(r13,r13)))))
    
    #--------------------------------
    #   BEST_SPOT
    #
    # Requires:
    #   xyz0    | [array] (3,) center of the sphere
    #   xyz2    | [array] (N,3) coordinates of the neighbor atoms
    #   cell    | [array] (3,3) cell matrix of shape. See MIC.
    #   r       | [float] radius of the sphere.
    #   trials  | [float] number of trials for the MC sampling.
    #
    # Returns:
    #   sp      | [array] (trials, 3) coordinates of all trial points
    #   sp_opt  | [array] (3,) coordinates of the best spot
    #   dist    | [float] distance to the clostest neighbor
    #   cone    | [float] cone angle in degrees.
    #
    # Notes:
    #   --
    #
    def best_spot(self,xyz0, xyz2, cell, r, trials):
        sp=[]
        for ii in range(trials):
            sp.append(self.sphere_point([xyz0[0], xyz0[1], xyz0[2]],r))
        PDM=self.pdm(np.array(sp), np.array(xyz2), cell)
        dist=np.amax([np.amin(ii) for ii in PDM])
        opt_idx=np.argmax([np.amin(ii) for ii in PDM])
        closest=xyz2[np.argmin(PDM[opt_idx])]
        cone=self.angle(xyz0, closest, sp[opt_idx],cell)
        return sp, sp[opt_idx], dist, cone
    
