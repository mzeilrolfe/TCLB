#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 30 22:42:20 2022

@author: mzeilrolfe
"""

import sys
sys.path.append('TCLB_tools/Python')
import CLB.VTIFile


fvti = 'TCLB/output/GA_4_VTK_P00_00151000.pvti'
vti = CLB.VTIFile.VTIFile(fvti, True)
vti = CLB.VTIFile.VTIFile3D(fvti, True)
phasefield = vti.get('PhaseField')

Boundary = vti.get('BOUNDARY')


xmax = Boundary.shape[0]-1
ymax = Boundary.shape[1]-1
zmax = Boundary.shape[2]-1
nw = 0.03

"""
Calculate Contact Area
"""
#velocity   = vti.get('U', vector=True)
count = 0 # track number of bounadry nodes next to liquid phase
count1 = 0
i=0
j=0
k=0
zmin = round(0.8/nw)
for k in list(range(zmin,zmax)):
     for j in list(range(ymax)):
         for i in list(range(xmax)):
             if Boundary[i,j,k] > 0.1:
                 if Boundary[i,j,k+1] == 0 and phasefield[i,j,k+1]>0.5:
                     count+=1
                     count1+=1
                 if Boundary[i+1,j,k] == 0 and phasefield[i+1,j,k]>0.5:
                     count+=1
                 if Boundary[i-1,j,k] == 0 and phasefield[i-1,j,k]>0.5:
                     count+=1
                 if Boundary[i,j+1,k] == 0 and phasefield[i,j+1,k]>0.5:
                     count+=1
                 if Boundary[i,j-1,k] == 0 and phasefield[i,j-1,k]>0.5:
                     count+=1        
An = nw**2 #Node Area (mm2)
print(round(count*An,4)) 
print(round(count1*An,4))
"""
Calculate Height of Centre of Mass
"""
i=0
j=0
k=0
klist=[]
for k in list(range(18,zmax)):
     for j in list(range(ymax)):
         for i in list(range(xmax)):
             if Boundary[i,j,k] < 1 and phasefield[i,j,k] > 0.5:
                 klist.append(k)
COMz2 = sum(klist)/len(klist)
print(round(COMz2*nw-1,3))
