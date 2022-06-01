#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 04:38:17 2022

@author: mzeilrolfe
"""


"""Due to the size of these vtk files, these had to be processed
1 by 1 and then the velocities were calculated and plotted in 
MATLAB
"""


import sys
sys.path.append('TCLB_tools/Python')
import CLB.VTIFile

 
fvti = 'TCLB/output/P2_Plain_VTK_P00_00181000.pvti'
vti = CLB.VTIFile.VTIFile(fvti, True)

vti = CLB.VTIFile.VTIFile3D(fvti, True)

phasefield = vti.get('PhaseField')

Boundary = vti.get('BOUNDARY')


nw = 0.05
"""
Calculate Droplet position down wedge
"""
i=0
j=0
k=0
xlist=[]
for k in list(range(19,59)):
    for j in list(range(415)):
        for i in list(range(99)):
            
            if Boundary[i,j,k]==0:
                if phasefield[i,j,k] > 0.5:
                    xlist.append(j)

BF = max(xlist)*nw #Bubble Front Positions
BB = min(xlist)*nw #Bubble Back Positions
print(BB)
print(BF)