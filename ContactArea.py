#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 04:07:58 2022

@author: mzeilrolfe
"""

import sys
sys.path.append('TCLB_tools/Python')
import CLB.VTIFile


fvti = 'TCLB/output/P2_C1_3_VTK_P00_00001000.pvti'
vti = CLB.VTIFile.VTIFile(fvti, True)
vti = CLB.VTIFile.VTIFile3D(fvti, True)


phasefield = vti.get('PhaseField')

Boundary = vti.get('BOUNDARY')


"""
Calculate Contact Area
"""
#velocity   = vti.get('U', vector=True)
count = 0 # track number of bounadry nodes next to liquid phase
count1 = 0
i=0
j=0
k=0
for k in list(range(18,28)):
     for j in list(range(139)):
         for i in list(range(139)):
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
An = 0.05**2 #Node Area (mm2)
print(round(count*An,2)) 
print(round(count1*An,2))
"""
Calculate Height of Centre of Mass
"""
i=0
j=0
k=0
ksum=0
PF_sum=0
klist=[]
for k in list(range(18,100)):
     for j in list(range(139)):
         for i in list(range(139)):
             if Boundary[i,j,k] < 1 and phasefield[i,j,k] > 0.5:
                 ksum += (k)
                 PF_sum += phasefield[i,j,k]
                 klist.append(k)
COMz = ksum/PF_sum
COMz2 = sum(klist)/len(klist)
print(round(COMz2*0.05-1,3))