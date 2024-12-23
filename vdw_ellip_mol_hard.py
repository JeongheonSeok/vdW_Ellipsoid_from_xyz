#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 16:05:07 2024

@author: JeongHeonSeok
"""

from __future__ import division
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sys
import numpy as np
from numpy import linalg
from random import random

VDWR = {"H":1.20,"He":1.40,
        "Li":1.82,"Be":1.53,"B":1.92,"C":1.70,"N":1.55,"O":1.52,"F":1.47,"Ne":1.54,
        "Na":2.27,"Mg":1.73,"Al":1.84,"Si":2.10,"P":1.80,"S":1.80,"Cl":1.75,"Ar":1.88,
        "K":2.75,"Ca":2.31,"Sc":2.11,"Ni":1.63,"Cu":1.40,"Zn":1.39,"Ga":1.87,"Ge":2.11,"As":1.85,"Se":1.90,"Br":1.85,"Kr":2.02,
        "Rb":3.03,"Sr":2.49,"Pd":1.63,"Ag":1.72,"Cd":1.58,"In":1.93,"Sn":2.17,"Sb":2.06,"Te":2.06,"I":1.98,"Xe":2.16,
        "Cs":3.43,"Ba":2.68,"Pt":1.75,"Au":1.66,"Hg":1.55,"Tl":1.96,"Pb":2.02,"Bi":2.07,"Po":1.97,"At":2.02,"Rn":2.20,
        "Fr":3.48,"Ra":2.83,"U":1.86}   # https://en.wikipedia.org/wiki/Van_der_Waals_radius
N_sph = 100

class EllipsoidTool:
    """Some stuff for playing with ellipsoids"""
    def __init__(self): pass
    
    def getMinVolEllipse(self, P=None, tolerance=0.01):
        """
        from ant-trullo/smFiSH_software/EllipsoidTool.py
        https://github.com/ant-trullo/smFiSH_software/blob/main/EllipsoidTool.py
        
        Find the minimum volume ellipsoid which holds all the points
        
        Based on work by Nima Moshtagh
        http://www.mathworks.com/matlabcentral/fileexchange/9542
        and also by looking at:
        http://cctbx.sourceforge.net/current/python/scitbx.math.minimum_covering_ellipsoid.html
        Which is based on the first reference anyway!
        
        Here, P is a numpy array of N dimensional points like this:
        P = [[x,y,z,...], <-- one point per line
             [x,y,z,...],
             [x,y,z,...]]
        
        Returns:
        (center, radii, rotation)
        
        """
        (N, d) = np.shape(P)    # N:점 개수 , d:공간 차원
        d = float(d)
    
        # Q will be our working array
        Q = np.vstack([np.copy(P.T), np.ones(N)]) 
        QT = Q.T
        """
        QT = [[x1,y1,z1,...,1],
              [x2,y2,z2,...,1],
              ...]
        """
        
        # initializations
        err = 1.0 + tolerance
        u = (1.0 / N) * np.ones(N)  # np.array([1/N,1/N,...,1/N], 개수 N개)

        # Khachiyan Algorithm
        while err > tolerance:
            V = np.dot(Q, np.dot(np.diag(u), QT))                 # (Q \dot Q^T)/N
            M = np.diag(np.dot(QT , np.dot(linalg.inv(V), Q)))    # M the diagonal vector of an NxN matrix
            j = np.argmax(M)
            maximum = M[j]
            step_size = (maximum - d - 1.0) / ((d + 1.0) * (maximum - 1.0))
            new_u = (1.0 - step_size) * u
            new_u[j] += step_size
            err = np.linalg.norm(new_u - u)
            u = new_u

        # center of the ellipse 
        center = np.dot(P.T, u)
    
        # the A matrix for the ellipse
        A = linalg.inv(
                       np.dot(P.T, np.dot(np.diag(u), P)) - 
                       np.array([[a * b for b in center] for a in center])
                       ) / d
                       
        # Get the values we'd like to return
        U, s, rotation = linalg.svd(A)
        radii = 1.0/np.sqrt(s)

        #print(center)
        #print(radii)
        #print(rotation)
        
        return (center, radii, rotation)

    def getEllipsoidVolume(self, radii):
        """Calculate the volume of the blob"""
        return 4./3.*np.pi*radii[0]*radii[1]*radii[2]

    def plotEllipsoid(self, center, radii, rotation, ax=None, plotAxes=False, cageColor='b', cageAlpha=0.2):
        """Plot an ellipsoid"""
        make_ax = ax == None
        if make_ax:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            
        u = np.linspace(0.0, 2.0 * np.pi, 100)
        v = np.linspace(0.0, np.pi, 100)
        
        # cartesian coordinates that correspond to the spherical angles:
        x = radii[0] * np.outer(np.cos(u), np.sin(v))
        y = radii[1] * np.outer(np.sin(u), np.sin(v))
        z = radii[2] * np.outer(np.ones_like(u), np.cos(v))
        # rotate accordingly
        for i in range(len(x)):
            for j in range(len(x)):
                [x[i,j],y[i,j],z[i,j]] = np.dot([x[i,j],y[i,j],z[i,j]], rotation) + center
    
        if plotAxes:
            # make some purdy axes
            axes = np.array([[radii[0],0.0,0.0],
                             [0.0,radii[1],0.0],
                             [0.0,0.0,radii[2]]])
            # rotate accordingly
            for i in range(len(axes)):
                axes[i] = np.dot(axes[i], rotation)
    
    
            # plot axes
            for p in axes:
                X3 = np.linspace(-p[0], p[0], 100) + center[0]
                Y3 = np.linspace(-p[1], p[1], 100) + center[1]
                Z3 = np.linspace(-p[2], p[2], 100) + center[2]
                ax.plot(X3, Y3, Z3, color=cageColor)
    
        # plot ellipsoid
        ax.plot_wireframe(x, y, z,  rstride=4, cstride=4, color=cageColor, alpha=cageAlpha)
        
        if make_ax:
            plt.show()
            plt.close(fig)
            del fig
            
def generate_points_on_sphere(x_c, y_c, z_c, atom, N=N_sph):
    # with golden spiral method
    r = VDWR[atom]
    points = []
    phi = (1 + np.sqrt(5)) / 2  # 황금비
    for i in range(N):
        z = 1 - (2 * (i + 0.5) / N)  # z축 좌표
        theta = 2 * np.pi * i / phi  # 경도
        x = np.sqrt(1 - z**2) * np.cos(theta)
        y = np.sqrt(1 - z**2) * np.sin(theta)
        # 구 좌표 -> 직교 좌표 변환
        points.append([x_c + r * x, y_c + r * y, z_c + r * z])
    return points
            
def read_xyz(xyzdir):
    atoms = []
    coords = []
    with open(xyzdir, 'r') as file:
        lines = file.readlines()
    i = 0
    for line in lines:
        i += 1
        if i<3:
            continue
        spline = line.split()
        atoms.append(spline[0])
        coords.append([float(spline[1]),float(spline[2]),float(spline[3])])
    
    return (atoms, coords)

def getvdwCoords(xyzdir,N=N_sph):
    (atoms, orgcoords) = read_xyz(xyzdir)
    vdwcoords = []
    for i in range(len(atoms)):
        points = generate_points_on_sphere(orgcoords[i][0],orgcoords[i][1],orgcoords[i][2],atoms[i],N)
        for point in points:
            vdwcoords.append(point)
    return vdwcoords

def getMinVolmolEllipse(xyzdir,tolerance=0.01,N=N_sph):
    (atoms, coords) = read_xyz(xyzdir)
    coords = np.array(coords)
    ET = EllipsoidTool()
    (center, radii, rotation) = ET.getMinVolEllipse(coords,tolerance)
    return (center, radii, rotation)

def getMinVolvdwEllipse(xyzdir,tolerance=0.01,N=N_sph):
    vdwcoords = np.array(getvdwCoords(xyzdir,N))
    ET = EllipsoidTool()
    (center, radii, rotation) = ET.getMinVolEllipse(vdwcoords,tolerance)
    return (center, radii, rotation)
    
if __name__ == "__main__":
    xyzdir = "C:/Users/JeongHeonSeok/Desktop/DIPEA.xyz"
    tolerance = 0.01
    N_sph = 100
    (nuccenter, nucradii, nucrotation) = getMinVolmolEllipse(xyzdir,tolerance = tolerance, N = N_sph)
    print(nucradii)
    (center, radii, rotation) = getMinVolvdwEllipse(xyzdir,tolerance = tolerance, N = N_sph)
    print(radii)
    print((radii[0]*radii[1]*radii[2])**(1/3))
    
    
"""
if __name__ == "__main__":
    # make 100 random points
    P = np.reshape([random()*100 for i in range(300)],(100,3))

    # find the ellipsoid
    ET = EllipsoidTool()
    (center, radii, rotation) = ET.getMinVolEllipse(P, .01)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # plot points
    ax.scatter(P[:,0], P[:,1], P[:,2], color='g', marker='*', s=100)

    # plot ellipsoid
    ET.plotEllipsoid(center, radii, rotation, ax=ax, plotAxes=True)
    
    plt.show()
    plt.close(fig)
    del fig
"""
