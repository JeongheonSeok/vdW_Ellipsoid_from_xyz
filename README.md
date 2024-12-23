# vdW_Ellipsoid_from_xyz
Code to find the minimum ellipsoid that encloses a molecule (where the volume is defined by the van der Waals radius) given in the form of a .xyz file.

## vdw_ellip_mol_hard.py
1. Read the coordinates of each atom that makes up the molecule from the ~.xyz file.
2. For each atom, evenly generate N points on the sphere of the van der Waals radius that surrounds it. ((Total number of atoms) * N points generated)
3. Find a minimum volume ellipsoid that encloses all points.

Since the sphere is approximated by N points, the larger the N is, the more accurate the result is obtained, but it takes longer.

It also used a Brute force algorithm that generates all the maximum points, not the minimum points required for actual ellipsoid calculations. So it still needs to be improved in terms of speed.

If anyone has improved it, I would greatly appreciate it if you could let me know.
