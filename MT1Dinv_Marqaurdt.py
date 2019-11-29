
#------------------------------------------------------------------------------
#                              MT1Dinv_Marqaurdt.py
# 
#  - This python program Produce the input data, initial model and run a 1D MT
#    Inversion Program written by Associate Professor Hideki Mizunaga (kyushu 
#    University,Japan). The MT 1D Inversion FORTRAN code was written using 
#    the modefied Marquardt method. 
# 
#                            HOW THE PROGRAM WORKS:
#                               ----------
#   - Put this code, the edifile and the Fortran executable (1dMarquardt.exe) in
#     the same directory. 
#   - Then, Frequency and Impedance components are extracted from the EDI
#     file.
#   - And the code calculates the apparent resistivity(TE,TM and DET modes) from the 
#     impedance components.
#   - The program then writes out the input data and initial model for the 1D inversion.
#     The default is the DET mode. The DET is calculated with the following formula
#     (Ranganayaki,1984; Vozoff, 1991)
      
#                  RHO_DET = 0.2*Period*|zxx*zyy-zxy*zyx|.
#                      
#    The executable (1dMarquardt.exe) is produced from the Fortran code(1DinvNLSQ.f90)
#    with the command line, cmd or cygwin. Using the following command 
#               gfortran -o 1dMarquardt.exe 1DinvNLSQ.f90
#
#   * The user should input 
#     - Initial model parametres
#     - target depth to be plotted 
#     - Fortran Excutable file 
#     - Location of Project Area (Not more than 40 charcters)
#     - Station name
#
#  << References>>
#   - Mizunaga, H., (1992) Ver. 1.01  , MT 1D Inversion code.
#   - Levenberg, K., (1944). A method for the solution of certain nonlinear 
#     Problems in least squares, Quarterly of Applied Mathematics, 2, 164-168.
#   - Marqaurdt, D.W. (1963). An algorithm for least squares estimation of non-
#     linear parameters. Journal of the Society of Industrial and Applied
#     Mathematics, 11, 431-441.
#   - Ranganayaki, R. P. (1984). An interpretative analysis of magnetotelluric 
#     data, Geophysics, 49, 1730–1748.
#   - Vozoff, K. (1991). The magnetotelluric method. In: M. N. Nabighian (ed),
#     Electromagnetic Methods in Applied Geophysics, volume 2 - Application,  
#     chapter 8, Society Exploration Geophysicists, 641-711.  
#
# By Biruk A. Cherkose @2017
#  email: 201990083@uaeu.ac.ae                                   -----
#                                                                      

import math
from subprocess import Popen, PIPE
import matplotlib.pyplot as plt 
import numpy as np
from sklearn.metrics import mean_squared_error # imported to get a RMS error function
from math import sqrt

#--------------INPUTs from the user--------------------------------------------------
station = 'AYR-301' # sation name
Location = 'Butajira, MER, ETHIOPIA' # Location of Project Area  
RhoInitial = [18,16,6,3,14,4] # [Rho1,Rho2,..Rho(n)] Initial model resistivity values
ThiInitial = [100,100,100,100,100] # [Thi1,Thi2,..Thi(n-1)] Initial model layers Thickness
exefilename = '1dMarquardt.exe' # Fortran Excutable file name
Targetdepth= 700 # Depth to plot in meters 


#---------------------------------------------------------------------------------
NoParameter= len(RhoInitial)+len(ThiInitial)
NumLayers = len(RhoInitial)
f = open("{}.edi".format(station))
EDIFile = f.readlines() 
invresult = "{}_1Dout.txt".format(station)# Output of Inversion

#-----------------Extracting No.Freq, lat and Long ----------------------------
nofreq1 = [] 
lat1=[] 
long1 = []
for index, line in enumerate(open("{}.edi".format(station))):
    if 'NFREQ=' in line:       
        nofreq1.append(line.split()[0])
    if 'REFLAT=' in line:
        lat1.append(line.split()[0])
    if 'REFLONG=' in line:
        long1.append(line.split()[0])
       
N = nofreq1[0]
Nofreq = N[6:]
nlat = lat1[0]         
#print nlat
LAT = nlat[7:]
nlong = long1[0]         
#print nlong
LONG = nlong[8:]
NoFrequency = Nofreq # No. of frequency from the EDI
#-------------------------Extracting Frequency---------------------------------      
for i, line in enumerate(EDIFile):
    if '>FREQ' in line:
        freqstart = i
    elif '>' in line:
        freqstop = i
        
freq1 = []

for j in range(freqstart+1, freqstop):
    raw_freq = [k for k in EDIFile[j].split('\n')]
    freq1.extend([l for l in raw_freq[0].split()]) 

ff = freq1[:int(Nofreq)] 

freq = []
for item in ff:
    freq.append(float(item))

#-----------------------ZXXReal------------------------------------------------
for i, line in enumerate(EDIFile):
    if '>ZXXR' in line:
        impstart = i
    elif '>' in line:
        impstop = i
#step = impstop-impstart
        
ZXXR1 = []

for j in range(impstart+1, impstop):
    raw_imp = [k for k in EDIFile[j].split('\n')]
    ZXXR1.extend([l for l in raw_imp[0].split()]) 

ZXXReal = ZXXR1[:int(Nofreq)]
#-----------------------ZXXImag------------------------------------------------
for i, line in enumerate(EDIFile):
    if '>ZXXI' in line:
        impstart = i
    elif '>' in line:
        impstop = i
#step = impstop-impstart
        
ZXXI1 = []
for j in range(impstart+1, impstop):
    raw_imp = [k for k in EDIFile[j].split('\n')]
    ZXXI1.extend([l for l in raw_imp[0].split()]) 

ZXXImag = ZXXI1[:int(Nofreq)]

#-----------------------ZXYReal------------------------------------------------
for i, line in enumerate(EDIFile):
    if '>ZXYR' in line:
        impstart = i
    elif '>' in line:
        impstop = i
#step = impstop-impstart
        
ZXYR1 = []

for j in range(impstart+1, impstop):
    raw_imp = [k for k in EDIFile[j].split('\n')]
    ZXYR1.extend([l for l in raw_imp[0].split()]) 

ZXYReal = ZXYR1[:int(Nofreq)]

#-----------------------ZXYImag------------------------------------------------
for i, line in enumerate(EDIFile):
    if '>ZXYI' in line:
        impstart = i
    elif '>' in line:
        impstop = i
#step = impstop-impstart
        
ZXYI1 = []

for j in range(impstart+1, impstop):
    raw_imp = [k for k in EDIFile[j].split('\n')]
    ZXYI1.extend([l for l in raw_imp[0].split()]) 

ZXYImag = ZXYI1[:int(Nofreq)]

#-----------------------ZYXReal------------------------------------------------
for i, line in enumerate(EDIFile):
    if '>ZYXR' in line:
        impstart = i
    elif '>' in line:
        impstop = i
#step = impstop-impstart
        
ZYXR1 = []

for j in range(impstart+1, impstop):
    raw_imp = [k for k in EDIFile[j].split('\n')]
    ZYXR1.extend([l for l in raw_imp[0].split()]) 

ZYXReal = ZYXR1[:int(Nofreq)]

#-----------------------ZYXImag------------------------------------------------
for i, line in enumerate(EDIFile):
    if '>ZYXI' in line:
        impstart = i
    elif '>' in line:
        impstop = i
#step = impstop-impstart
        
ZYXI1 = []

for j in range(impstart+1, impstop):
    raw_imp = [k for k in EDIFile[j].split('\n')]
    ZYXI1.extend([l for l in raw_imp[0].split()]) 

ZYXImag = ZYXI1[:int(Nofreq)]

#-----------------------ZYYReal------------------------------------------------
for i, line in enumerate(EDIFile):
    if '>ZYYR' in line:
        impstart = i
    elif '>' in line:
        impstop = i
#step = impstop-impstart
        
ZYYR1 = []

for j in range(impstart+1, impstop):
    raw_imp = [k for k in EDIFile[j].split('\n')]
    ZYYR1.extend([l for l in raw_imp[0].split()]) 

ZYYReal = ZYYR1[:int(Nofreq)]


#-----------------------ZYYImag------------------------------------------------
for i, line in enumerate(EDIFile):
    if '>ZYYI' in line:
        impstart = i
    elif '>' in line:
        impstop = i
#step = impstop-impstart
        
ZYYI1 = []

for j in range(impstart+1, impstop):
    raw_imp = [k for k in EDIFile[j].split('\n')]
    ZYYI1.extend([l for l in raw_imp[0].split()]) 

ZYYImag = ZYYI1[:int(Nofreq)]

#----------------------------------------------------------------------------
ZXXR = []
for item in ZXXReal:
    ZXXR.append(float(item))
ZXXI = []
for item in ZXXImag:
    ZXXI.append(float(item))

ZXYR = []
for item in ZXYReal:
    ZXYR.append(float(item))
ZXYI = []
for item in ZXYImag:
    ZXYI.append(float(item))

ZYXR = []
for item in ZYXReal:
    ZYXR.append(float(item))
ZYXI = []
for item in ZYXImag:
    ZYXI.append(float(item))

ZYYR = []
for item in ZYYReal:
    ZYYR.append(float(item))
ZYYI = []
for item in ZYYImag:
    ZYYI.append(float(item))

#------------------------------------------------------------------------------
zxx = []; zyy = []; zxy = []; zyx = []
ZXXc = []; ZYYc = []; ZXYc = []; ZYXc = []

period = []

RHOxy = []; RHOyx = []; RHOxx = []; RHOyy = []
phyYX=[]; phyXY = []; phyXX= []; phyYY= []

Rho_DET = []

for i in range(len(freq)):
    ZXXc.append(complex(ZXXR[i],ZXXI[i]))
    zxx.append(abs(ZXXc[i]))
    ZXYc.append(complex(ZXYR[i],ZXYI[i]))
    zxy.append(abs(ZXYc[i]))
    ZYXc.append(complex(ZYXR[i],ZYXI[i]))
    zyx.append(abs(ZYXc[i]))
    ZYYc.append(complex(ZYYR[i],ZYYI[i]))
    zyy.append(abs(ZYYc[i]))
    
    phyYX.append(math.degrees(math.atan(ZYXc[i].imag/ZYXc[i].real)))
    phyXY.append(math.degrees(math.atan(ZXYc[i].imag/ZXYc[i].real)))
    phyXX.append(math.degrees(math.atan(ZXXc[i].imag/ZXXc[i].real)))
    phyYY.append(math.degrees(math.atan(ZYYc[i].imag/ZYYc[i].real)))
    
    period.append(1/freq[i])
    
    RHOxy.append(0.2*period[i]*abs(ZXYc[i]**2))
    RHOyx.append(0.2*period[i]*abs(ZYXc[i]**2))
    RHOxx.append(0.2*period[i]*abs(ZXXc[i]**2))
    RHOyy.append(0.2*period[i]*abs(ZYYc[i]**2))
    Rho_DET.append((0.2*abs(ZXXc[i]*ZYYc[i]-ZXYc[i]*ZYXc[i]))/freq[i])
#print  Rho_DET
#print freq
    
LOC = "LOC {}\n".format(Location)
NoFreq = "{}\n".format(NoFrequency)
NoParam = "{}\n".format(NoParameter)
#print NoParam
#------------Arrange frequency from small to large value----------------------
frequency = list(reversed(freq))
Determinant = list(reversed(Rho_DET))
  
Inputfile = '{}_1Dinput.hm'.format(station)    
with open(Inputfile.format(station), 'w') as text_file:
    text_file.write(LOC)
    text_file.write(NoFreq) # number of frequencies
    for j in range(len(freq)):
            text_file.write('{:7.7s}\t{:7.7s}\t1\n'.format(str(frequency[j]),
                         str(Determinant[j])))
    text_file.write(NoParam) # number of layer
    for x in range(len(RhoInitial)):
        text_file.write(('{:7.7s}\t0\n'.format(str(RhoInitial[x]))))
    for x in range(len(ThiInitial)):
        text_file.write(('{:7.7s}\t0\n'.format(str(ThiInitial[x]))))   
    

#----------RUN Fortan exe file from Python IDE or Editor-----------------------

LOC = "LOC {}\n".format(Location)
input = open(Inputfile, 'rb').read()

running_proc = Popen([exefilename], stdout=PIPE, stdin=PIPE, stderr=PIPE)
out, err = running_proc.communicate(input=input)
outfile = open(invresult, "w")
outfile.write(out.decode())
outfile.close()

##---------------------Part II: Plotting Inversion Result-----------------------

##---------------------Extracting from inv. output file-------------------------
freq2 = []
Obs1 = []
Cal1 = []
Par1 = [] 
Par2 = [] 
N = 8 + (2*NumLayers-1)   
for index, line in enumerate(open(invresult)):
    if index > 6 and index <= 6+(2*NumLayers-1):
        Par2.append(line.split()[2])
        Par1.append(line.split()[3])
    elif index > N:
        freq2.append(line.split()[1])
        Obs1.append(line.split()[2])
        Cal1.append(line.split()[3])        
    else:
        continue


#-------------------Converting values from string to float----------------------
fff = []
Obs = []
Cal = []
Par = []
Pari = []

for item in freq2:
    fff.append(float(item))
for item in Obs1:
    Obs.append(float(item))
for item in Cal1:
    Cal.append(float(item))
for item in Par1:
    Par.append(float(item))
for item in Par2:
    Pari.append(float(item))
##-------------------------------------------------------------------------------
RHOinitial = []
THIinital= []
RHOfinal = []
THIfinal = []
RHOfinal = Par[:len(Par)/2+1] # Separate Rho from Thickness in the Parameter(par)
THIfinal = Par[len(Par)/2+1:] # Separate Thickness from Rho in the Parameter(par)
THIinitial = Pari[len(Pari)/2+1:]
RHOinitial = Pari[:len(Pari)/2+1]
Resi =np.repeat(RHOinitial,2)
Resf =np.repeat(RHOfinal,2)


y1 = [THIfinal[:i+1] for i in range(len(THIfinal))] # Group list to add layers 
y2 = [reduce(lambda x, y: x+y,l) for l in y1] # Add Layers
y3 = [x+.1 for x in y2]
Depth = sorted([0,Targetdepth] + y2 + y3) # Calculates depth for ploting

y1i = [THIinitial[:i+1] for i in range(len(THIinitial))] # Group list to add layers 
y2i = [reduce(lambda x, y: x+y,l) for l in y1i] # Add Layers
y3i = [x+.1 for x in y2i]
Depthf = sorted([0,Targetdepth] + y2 + y3) # Calculates depth for ploting final model
Depthi = sorted([0,Targetdepth] + y2i + y3i) # Calculates depth for ploting initial model
RMSerr=sqrt(mean_squared_error(Obs,Cal)) 

###------------------Capturing the No.of Iteration and RMS------------------------
iterno = []
for index, line in enumerate(open(invresult)):
    if index == 2:
       iterno.append(line.split()[5])
Iterno = []
IterationNo = Iterno
for item in iterno:
    Iterno.append(int(item))
#print Iterno
#print "%.2f"%RMSerr
print RHOfinal
print THIfinal
##------------------ Calculating Period from frequency----------------------------
period = []
P = []
observed = []
calculated = []
for i in range(len(fff)):
    period.append(1/fff[i])

#----------------------Ordering the period from Small to large Value -----------
P = list(reversed(period))  
observed = list(reversed(Obs))
calculated = list(reversed(Cal))
#-----------------------------------------------------------------------------
MisfitTitleID = "Observed vs. Calculated"
ModelTitleID = "Model-{}\nIteration = {} and RMS = {}".format(station,Iterno[0],"%.2f"%RMSerr) 
Title = MisfitTitleID.format(station) # Title for Mifit curve to be ploted
Modelname = ModelTitleID.format(station) # Title for 1-D Model to be ploted
Fig1 = "Model{}.jpg"

Misfit = Fig1.format(station) 

x = range(1,len(RHOfinal)+1) # X is used in the making of the model tabel


#--------------------------PLOT-----------------------------------------------
plt.style.use('seaborn-white')
plt.figure(figsize=(17,8))
#print(plt.style.available)
#-------------------------Plot Model-------------------------------------------
plt.subplot(1,2,1)

#plt.semilogx(Resi,Depthi,color = 'gray',label='Initial Model',lw=1,
#         linestyle = '--')
plt.semilogx(Resf,Depthf,color = 'gray',label='Final Model',lw=4)
plt.xscale('log')
#plt.yscale('log')
plt.gca().invert_yaxis()
plt.xticks(size=15,weight = 'bold')
plt.yticks(size=15,weight = 'bold')
plt.legend(fontsize=15)
plt.xlim(min(Resf*.15),1000)
plt.grid(True,which ='both')
plt.ylim(Targetdepth,1)
plt.title(Modelname, size=18,weight = 'bold')
plt.xlabel(s='Resistivity [ohm.m]',size=18,weight = 'bold')
plt.ylabel(s='Depth [m]',size=18,weight = 'bold')

#---------------Plot Misfit between Observed and Calculated--------------------

plt.subplot(1,2,2)

plt.plot(P,observed, label='Observed',color = 'k', marker = '+',linestyle = ' ')
plt.plot(P,calculated,label='Calculated',color = 'gray')

plt.xscale('log')
plt.yscale('log')
plt.xticks(size=15,weight = 'bold')
plt.yticks(size=15,weight = 'bold')
minval=min( min(observed,calculated))
maxval=max(max(observed,calculated))
plt.ylim([minval/10,maxval*10])
plt.xlim(0.5*min(period),2*max(period))
plt.xlabel(s='Period [s]',size=18,weight = 'bold')
plt.ylabel(s='App. Res. [ohm.m]',size=18,weight = 'bold')
plt.grid(which ='minor',axis = 'both',alpha=.5)
plt.title(Title,size=18,weight = 'bold')
plt.legend(fontsize=15, loc = 'best')
plt.savefig(Misfit,dpi = 300)

#'''
#plotModel tabel
#bbox is used for positioning the tabel
#'''
#plt.subplot(1,3,3)
#Layerno=['Layer Number'] + x
#Resfinaltab = ['RES'] + RHOfinal
#Thitab = ['THI'] + THIfinal + [' ']
#table_vals = [Layerno,Resfinaltab,Thitab]
#Table = plt.table(cellText = table_vals,
#                loc =  'best',
#                 bbox =[0,0.0,1,1])
#
#plt.xlim(0,2.5)
#plt.ylim(0,0.5)
#plt.xticks(size=0)
#plt.yticks(size=0)

plt.show()
