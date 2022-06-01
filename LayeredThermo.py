# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import sys
from vtk import *
from numpy import zeros
import matplotlib
import math
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
import matplotlib.pyplot as plt
plt.rcParams.update(plt.rcParamsDefault)
import pandas as pd
import numpy as np
import matplotlib
from matplotlib import rc
sys.path.append('TCLB_tools/Python')
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('text', usetex=True)
plt.close('all')

def fetchData( myFile ):
    a = vtkXMLPImageDataReader()
    a.SetFileName( myFile )
    a.Update()
    #f = vtkCellDataToPointData(f)
    #f.PointData.values

    myDomain = a.GetOutput().GetBounds()
    nx = myDomain[1]
    ny = myDomain[3]
    nz = myDomain[5]

    #f = servermanager.Fetch( f )
    f = a.GetOutput().GetCellData()
    return a, f, int(nx), int(ny), int(nz)
#_____________________________________________________________________________________________________
#List the file VTK and LOG files that you want to analyse and then we set up the arrays to store the information.

myFile    = r'TCLB/output/LayeredThermocapillaryFlow_v2_majidiPaper_K5_VTK_P00_00400000.pvti'   # sys.argv[1]
myDatFile = r'TCLB/output/LayeredThermocapillaryFlow_v2_majidiPaper_K5_Log_P00_00100000.csv'



myImg, myData, nx, ny, nz = fetchData( myFile )
myDat = pd.read_csv( myDatFile )

PF = myData.GetArray('PhaseField')
U  = myData.GetArray('U')
Temp  = myData.GetArray('T')
num = PF.GetNumberOfTuples()

pf = zeros([nx,ny,nz])
ux = zeros([nx,ny,nz])
uy = zeros([nx,ny,nz])
T  = zeros([nx,ny,nz])
myDatX = np.arange(nx)-myDat.myL[0] 
myDatY = np.arange(ny)-myDat.MIDPOINT[0]
XX, YY = np.meshgrid(myDatX, myDatY)


#_____________________________________________________________________________________________________
#Import The Data
for n in range(num):
    ind = myImg.GetOutput().GetCell(n).GetBounds()
    x = int(ind[0])
    y = int(ind[2])
    z = int(ind[4])
    pf[x,y,z] = PF.GetTuple(n)[0]
    utmp = U.GetTuple(n)
    ux[x,y,z] = utmp[0]
    uy[x,y,z] = utmp[1]
    T[x,y,z] = Temp.GetTuple(n)[0]



#_____________________________________________________________________________________________________
#Generate Analytic solution
from math import sinh, cosh, cos, sin, pi
myDat.myL = myDat.myL
myDat.MIDPOINT = myDat.MIDPOINT - 1.0
k_top = myDat.k_h[0]
k_bot = myDat.k_l[0]
kstar = k_top/k_bot
T_0 = myDat.T_0[0]
sigma_T = myDat.sigma_T[0]
T_c = myDat.T_c[0]
T_h = myDat.T_h[0]

w = pi/myDat.myL[0]
a = (myDat.MIDPOINT[0]) * w
b = (myDat.MIDPOINT[0]) * w

f = 1.0 / ((kstar)*sinh(b)*cosh(a) + sinh(a)*cosh(b))
g = sinh(a)*f

h = (sinh(a)**2-a**2)*(sinh(b)**2-b**2) / \
    ( (sinh(b)**2-b**2)*(sinh(2.0*a)-2.0*a)  \
        + (sinh(a)**2-a**2)*(sinh(2.0*b)-2.0*b))

Ca1 = sinh(a) ** 2 / (sinh(a)**2 - a**2)
Ca2 = -1.0 * (myDat.MIDPOINT[0]) * a / (sinh(a)**2 - a**2)
Ca3 = (2*a - sinh(2*a) ) / (2.0*(sinh(a)**2 - a**2))      

Cb1 = sinh(b) * sinh(b) / (sinh(b)*sinh(b) - b*b)
Cb2 = -1.0 * (myDat.MIDPOINT[0]) * b / (sinh(b)*sinh(b) - b*b)
Cb3 = (-2*b + sinh(2*b) ) / (2.0*(sinh(b)*sinh(b) - b*b))

umax = -1.0 * (T_0 * sigma_T/myDat.Viscosity_l[0]) * g * h
jj = 0
xx = np.linspace(-myDat.myL[0]-0.5,myDat.myL[0]-0.5,nx)
yy = np.linspace(-(myDat.MIDPOINT[0]),(myDat.MIDPOINT[0]),ny)
UX_a = zeros([len(xx),len(yy)])
UY_a = zeros([len(xx),len(yy)])
T_a = zeros([len(xx),len(yy)])
for y in yy:
    ii = 0
    for x in xx:
        if y > 0:
            T_a[ii,jj] = ((T_c - T_h)* y + (kstar)*T_c*(myDat.MIDPOINT[0]-1)+T_h*(myDat.MIDPOINT[0]-1)) / ((myDat.MIDPOINT[0]-1) + (kstar)*(myDat.MIDPOINT[0]-1)) + \
                        T_0*f * sinh(a - y*w) * cos(w*x)
            T_a[ii,jj] = ((T_c - T_h)* y + (kstar)*T_c*(myDat.MIDPOINT[0])+T_h*(myDat.MIDPOINT[0])) / ((myDat.MIDPOINT[0]) + (kstar)*(myDat.MIDPOINT[0])) + \
                        T_0*f * sinh(a - y*w) * cos(w*x)
            UX_a[ii,jj]  = umax* ( (Ca1+(w)*(Ca2+Ca3*y))*cosh(w*y) +  (Ca3+(w)*Ca1*y)*sinh(w*y) ) *sin(w*x)
            
            UY_a[ii,jj]  = -(w)*umax*(Ca1*y*cosh(w*y) + (Ca2+Ca3*y)*sinh(w*y))*cos(w*x)
            
        elif y <= 0:
            T_a[ii,jj] = ((kstar)*(T_c-T_h)*y + (kstar)*T_c*(myDat.MIDPOINT[0]-1) + T_h*(myDat.MIDPOINT[0]-1)) / ((myDat.MIDPOINT[0]-1) + (kstar)*(myDat.MIDPOINT[0]-1)) + T_0*f*(sinh(a)*cosh(w*y) - (kstar)*sinh(w*y)*cosh(a) )*cos(w*x)
            T_a[ii,jj] = ((kstar)*(T_c-T_h)*y + (kstar)*T_c*(myDat.MIDPOINT[0]) + T_h*(myDat.MIDPOINT[0])) / ((myDat.MIDPOINT[0]) + (kstar)*(myDat.MIDPOINT[0])) + T_0*f*(sinh(a)*cosh(w*y) - (kstar)*sinh(w*y)*cosh(a) )*cos(w*x)
        
            UX_a[ii,jj]  = umax* ( (Cb1+(w)*(Cb2+Cb3*y))*cosh(w*y)+(Cb3+(w)*Cb1*y)*sinh(w*y)) *sin(w*x)
            
            UY_a[ii,jj]  = -(w)*umax*(Cb1*y*cosh(w*y) + (Cb2+Cb3*y)*sinh(w*y))*cos(w*x)
        ii += 1            
    jj += 1
x,y = np.meshgrid(xx, yy)





excludeN = 10




# ANALYTIC SOLUTION
from scipy import integrate
intx = integrate.cumtrapz(UY_a.T, x, axis=1, initial=0)[0]
inty = integrate.cumtrapz(UX_a.T, y, axis=0, initial=0)
psi1=intx-inty
intx=integrate.cumtrapz(UY_a.T,x,axis=1,initial=0)
inty=integrate.cumtrapz(UX_a.T,y,axis=0,initial=0)[:,0][:,None]
psi2=intx-inty

psi=0.5*(psi1+psi2)
rangeC = [np.amin(psi), np.amax(psi)]
levels = np.linspace(rangeC[0],rangeC[1],20)
#Analytical


#Simulation results
intx = integrate.cumtrapz(uy[:,:,int(nz/2)].T, XX, axis=1, initial=0)[0]
inty = integrate.cumtrapz(ux[:,:,int(nz/2)].T, YY, axis=0, initial=0)
psi1=intx-inty
intx=integrate.cumtrapz(uy[:,:,int(nz/2)].T,XX,axis=1,initial=0)
inty=integrate.cumtrapz(ux[:,:,int(nz/2)].T,YY,axis=0,initial=0)[:,0][:,None]
psi2=intx-inty

psi=0.5*(psi1+psi2)
rangeC = [np.amin(psi), np.amax(psi)]
print(rangeC)

"""
L2temp_num = 0.0
L2temp_den = 0.0
L2u_num = 0.0
L2u_den = 0.0
from math import sqrt
for x in range(nx):
    for y in range(ny):
        u = sqrt(ux[x,y,int(nz/2)]**2 + uy[x,y,int(nz/2)]**2)
        U = sqrt(UX_a[x,y]**2 + UY_a[x,y]**2)
        L2temp_num +=  ( (T[x,y,int(nz/2)] - T_a[x,y])**2 )
        L2temp_den +=  ( (T_a[x,y])**2 )
        L2u_num +=  ( (u - U)**2 )
        L2u_den +=  ( (U)**2 )
print(" L2_T = %.8lf, L2_u = %.8lf" %(sqrt(L2temp_num/L2temp_den), sqrt(L2u_num/L2u_den)))
"""
np.savetxt('dataTk5.csv',T_a,delimiter=',')
np.savetxt('datauxk5.csv',UY_a[::excludeN,::excludeN],delimiter=',')
np.savetxt('datauyk5.csv',UX_a[::excludeN,::excludeN],delimiter=',')

np.savetxt('resultTk5.csv',T[:,:,int(nz/2)],delimiter=',')
np.savetxt('resultuxk5.csv',ux[::excludeN,::excludeN,int(nz/2)],delimiter=',')
np.savetxt('resultuyk5.csv',uy[::excludeN,::excludeN,int(nz/2)],delimiter=',')