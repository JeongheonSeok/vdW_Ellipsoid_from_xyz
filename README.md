# vdW_Ellipsoid_from_xyz
Code to find the minimum ellipsoid that encloses a molecule (where the volume is defined by the van der Waals radius) given in the form of a .xyz file.

## vdw_ellip_mol_hard.py
1. Read the coordinates of each atom that makes up the molecule from the ~.xyz file.
2. For each atom, evenly generate N_sph points on the sphere of the van der Waals radius that surrounds it. ((Total number of atoms) * N_sph points generated)
3. Find a minimum volume ellipsoid that encloses all points.

Since the sphere is approximated by N_sph points, the larger the N_sph is, the more accurate the result is obtained, but it needs more time.

It also used a Brute force algorithm that generates all the maximum points, not the minimum points required for actual ellipsoid calculations. So it still needs to be improved in terms of speed.

The time complexity of the current algorithm is O(N_atom * N_sph) for the number of atoms N_atom and the number of spherical approximative points N_sph. Using an algorithm that approximates spherical points only for areas close to the actual elliptic surface, it is expected that O((N_atom * N_sph)^(2/3)) for sufficiently large molecules.

If anyone has improved it, I would greatly appreciate it if you could let me know.

The key algorithm for finding the ellipsoid is from ant-trullo/smFiSH_software/EllipsoidTool.py
  https://github.com/ant-trullo/smFiSH_software/blob/main/EllipsoidTool.py
