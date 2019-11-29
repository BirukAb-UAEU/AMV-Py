
# -*- coding: utf-8 -*-
#""" 
#    *** Input_model_ModEM.py ***
#             
#    This python (2.7) code produces an input model file (3Dmode.rho)for ModEM program.
#    ModEM is a 3D MT inversion code(Egbert and Kelbert, 2012; Meqbel, 2009; Kelbert et
#    al., 2014). 
#
#    Biruk A. Cherkose
#    Email: 201990083@uaeu.ac.ae
#    @2019
#"""

import math
from  decimal import Decimal


#--------------USER Input -------------------------------------
cell_dimenstion_X = 500
cell_dimenstion_Y = 500
No_cells_x = 34 # No. of cells in X direction (North)
No_cells_y = 36 # No. of cells in Y direction (East)
x_pad = 7 # No. of cells added in X direction (North) on one side
y_pad = 7 # No. of cells added in Y direction (East) on one side
Nz = 55  # No of layers in the vertical direction
XandY_pad_factor = 1.3
Z_inc_factor = 1.1
Z_thic = 20
Initial_Rho = 15
#-------------------------------------------------------------

LOGE_Rho = math.log1p(Initial_Rho)

Nx = No_cells_x + 2*x_pad
Ny = No_cells_y + 2*y_pad

xx = []
DIM_x = cell_dimenstion_X

while (DIM_x < (cell_dimenstion_X * (XandY_pad_factor**x_pad)-1)):
    DIM_x = DIM_x * XandY_pad_factor
    xx.append("%.3f"%DIM_x)


yy = []
DIM_y = cell_dimenstion_Y

while (DIM_y < (cell_dimenstion_Y * (XandY_pad_factor**y_pad)-1)):
    DIM_y = DIM_y * XandY_pad_factor
    yy.append("%.3f"%DIM_y)
#print yy

zz = []
Zthickness = Z_thic

while (Zthickness < (Z_thic*Z_inc_factor**(Nz-1)+1)):
    Zthickness = Zthickness * Z_inc_factor
    zz.append("%.3f"%Zthickness)

print len(zz)
xcells = []
xcells.append(str("%.3f"%cell_dimenstion_X))
ycells = []
ycells.append(str("%.3f"%cell_dimenstion_Y))
    
xxx = list(reversed(xx)) + (No_cells_x*xcells) + xx
yyy = list(reversed(yy)) + (No_cells_y*ycells) + yy
rho_line = []
rho_line.append(LOGE_Rho)

rho_X = Nx*rho_line

m1 = sum(Decimal(i) for i in xxx)/2
m2 = sum(Decimal(i) for i in yyy)/2
print m1
print m2
with open("3Dmodel.rho", 'w') as wsfile:
    wsfile.write(" # 3D MT model written by ModEM in WS format\n")
    wsfile.write(" {}  {}  {}  0 LOGE\n".format(Nx,Ny,Nz))
    for ix in range(len(xxx)):
        wsfile.write(" {} ".format(str(xxx[ix])))
    wsfile.write("\n")
    for iy in range(len(yyy)):
        wsfile.write(" {} ".format(str(yyy[iy])))
    wsfile.write("\n")
#    wsfile.write("  {}\t  ".format("%.3f"%Z_thic))
    for ii in range(len(zz)):
        wsfile.write("   {} ".format(str(zz[ii])))
    wsfile.write("\n")
    wsfile.write("\n")
    j = 1
    while (j < Nz+1):
        
        jj = 1
        while (jj <  Ny+1):         
            for irho in range(len(rho_X)):
                wsfile.write("{} ".format(str('%.5E'%Decimal("%.5f"%rho_X[irho]))))
            jj = jj + 1
            wsfile.write("\n")
        j = j + 1
        wsfile.write("\n")
    wsfile.write("-{}    -{}    0.000\n".format(m1,m2))
#    wsfile.write("\n")
    wsfile.write("0.0")
    
print len(xxx)
print len(yyy)
         

#text_file.write('{:7.7s}\t +{:7.7s}\t {:7.7s}\n'.format(str(period[j]),
#                             str(freq[j]),str(DET[j])))
        
