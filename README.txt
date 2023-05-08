_________________________________________________________

	CONE ANALYSIS TOOL
_________________________________________________________

Version: 05-08-23

This is a preliminary version of my cone analysis tool. It provides cone angles
for any desired coordination polyhedra within an adjustable cut-off. It is
fully PBC-corrected and applicable to any triclinic lattices. Projection from
crystallographic -> orthogonal basis has to be done manually by now (see 'mic').
Comes with an example box of amorphous silicon nitride (a-Si3N4).
Detailed documentation is provided below.

Requirements:
- numpy
- random

Functions:
- mic
- pdm
- nl
- read_xyz
- sphere_point
- angle
- best_spot

--------------------------------
   MINIMUM IMAGE VECTOR (MIC)

 Requires:
   xi      | [array], shape (3,) coordinates of atom i
   xj      | [array], shape (3,) coordinates of atom j
   cell    | [array], shape (3,3) cell matrix of shape
           |   cell=[[ x,  0,  0],
                     [ xy, y,  0],
                     [ xz, yz, z]]

 Returns:
   minimum image distance, minimum image vector  

 Notes:
   Cell follows the lammps semantics, see 
   https://docs.lammps.org/Howto_triclinic.html
   
   Treatment of triclinic boxes is possible by projecting the box in an 
   orthogonal basis with the following requirements:
   
    1. All the basis vectors start in the origin.
    2. a is parallel to the x-axis and points in positive direction.
    3. b lies in the xy-plane with a positive y-component.
    4. c has a positive z-component.
    5. a,b,c form a right-handed basis. 
    
   This leads to 
      C=[[ x,  0,  0],
         [ xy, y,  0],
         [ xz, yz, z]]
         
    Cell parameters a, b, c, alpha, beta, gamma are projected from a 
    (non-orthognal) crystallographic basis into the required basis via
    
        x   = a
        xy  = b*cos(gamma)
        y   = b*sin(gamma)
        xz  = c*cos(beta)
        yz  = xy*c*cos(alpha)
        z   = sqrt(c**2-xz**2-yz**2)   


--------------------------------
   PAIR DISTANCE MATRIX (PDM)

 Requires:
   xyz1 | array/list (3,N) coordinates of the particles of interest
   xyz2 | array/list (3,M) coordinates of the reference particles.

 Returns:
   pair distance matrix (N,M)  

 Notes:
   Uses two distinct coord vectors. Computes the neighborlist for the atoms
   in xyz1 with respect to the ones in xyz2. Allows to get the PDM i.e. for 
   subsets of an atomistic system.
   If xyz1==xyz2, the computation is accelerated by computing only the upper
   half of the pdm.   

--------------------------------
   NEIGHBOR LIST (NL)

 Requires:
   xyz1      | array/list (3,N) coordinates of the particles of interest
   xyz2      | array/list (3,M) coordinates of the target particles.
   r_cut     | float () cut-off radius for 1st coordination shell
   cell      | [array], shape (3,3) cell matrix of shape. See MIC.
   R_ABS_MAT | array/list (N,M) if already given: The PDM. (optional)
 
 Returns:
   atom indices, neighbor indices (both zero bound)

 Notes:
   Uses two distinct coord vectors. Computes the neighborlist for the atoms
   in xyz1 with respect to the ones in xyz2. Allows to get the NL i.e. for 
   subsets of an atomistic system.

   Can use external PDM.

--------------------------------
   READ XYZ

 Requires:
   filename | (string) self explaining
 
 Returns:
   x, y, z
 Notes:
   only single frame, no trajectory

--------------------------------
   SPHERE POINT

 Requires:
   xyz0 | [array] (3,) center of the sphere
   r    | [float]  radius of the sphere
 
 Returns:
   [x, y, z] of the sphere point
 Notes:
   --

--------------------------------
   ANGLE

 Requires:
   xyz1    | [array] (3,) coordinates of the center atom
   xyz2    | [array] (3,) coordinates of the neighbor 1
   xyz3    | [array] (3,) coordinates of the neighbor 2
   cell    | [array] (3,3) cell matrix of shape. See MIC.

 Returns:
   angle (float) in degrees

 Notes:
   --

--------------------------------
   BEST_SPOT

 Requires:
   xyz0    | [array] (3,) center of the sphere
   xyz2    | [array] (N,3) coordinates of the neighbor atoms
   cell    | [array] (3,3) cell matrix of shape. See MIC.
   r       | [float] radius of the sphere.
   trials  | [float] number of trials for the MC sampling.

 Returns:
   sp      | [array] (trials, 3) coordinates of all trial points
   sp_opt  | [array] (3,) coordinates of the best spot
   dist    | [float] distance to the clostest neighbor
   cone    | [float] cone angle in degrees.

 Notes:
   --
