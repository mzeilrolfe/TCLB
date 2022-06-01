#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 04:07:58 2022

@author: mzeilrolfe
"""

import sys
sys.path.append('TCLB_tools/Python')
import CLB.VTIFile

fvti = 'TCLB/output/2D_noT_VTK_P00_04000200.pvti'
vti = CLB.VTIFile.VTIFile(fvti, True)
vti = CLB.VTIFile.VTIFile3D(fvti, True)


phasefield = vti.get('PhaseField')

Boundary = vti.get('BOUNDARY')

myshape = Boundary.shape
"""
Calculate Contact Area
"""
#velocity   = vti.get('U', vector=True)
count = 0 # track number of bounadry nodes next to liquid phase
pflim = 0.5
count1 = 0
i=3
j=0
k=0
for k in list(range(99)):
     for j in list(range(159)):
             
             if Boundary[i,j,k] > 0.1:
                 if Boundary[i,j,k+1] == 0 and phasefield[i,j,k+1]>pflim:
                     count+=1
ln = 0.05 #Node Area (mm2)
print(round(count*ln,4)) 
"""
Calculate Height of Centre of Mass
"""
i=3
j=0
k=0
ksum=0
PF_sum=0
klist=[]
for k in list(range(99)):
     for j in list(range(159)):
             if Boundary[i,j,k] < 1 and phasefield[i,j,k] > pflim:
                 ksum += (k)
                 PF_sum += phasefield[i,j,k]
                 klist.append(k)
COMz = ksum/PF_sum
COMz2 = sum(klist)/len(klist)
print(round(COMz2*0.05-1,4))
print(len(klist))
