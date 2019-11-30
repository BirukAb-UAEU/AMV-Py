"""
                Ediz_01.py
                ------
  - This Python (v2.7) program plots Apparent resistivity and Phase curves
  - Calculates and plots the Swift-Skew (κ) and Phase-Sensitive skew (η) 
  - file (EDI format == imp(z)).
  - Produce *.mt file that conatins caculated Apparent resistivity and Phase values
    and Swift-Skew (κ) and Phase-Sensitive skew (η).
    
HOW THE PROGRAM WORKS
   - Put the path of edi file (name_of_folder) in the USER INPUT section
   - Input the sation name (EDI file name)
   - Then RUN
   
OUTPUT after Runing the Program
   - JPG files for MT Curves and Swift-Skew (κ) and Phase-Sensitive skew (η)
   - *.mt file
Author: Biruk A. Cherkose, @ 2019, Email: 201990083@uaeu.ac.ae
"""

#---USER INPUT Section-------------
station = 'A01'
name_of_folder = r'.' 
#------------------------------------
#=======================================================


import os
import csv
import math
import matplotlib.pyplot as plt

#---INPUT-EDI file name and Folder---
station = 'A01'
name_of_folder = r'.' 
#------------------------------------

f = open("{}.edi".format(station))
EDIFile = f.readlines() 

#------ NFREQ, LAT and LONG ------------

nofreq1 = []; lat1=[]; long1 = [] 

for index, line in enumerate(open("{}.edi".format(station))):
    if 'NFREQ=' in line:       
        nofreq1.append(line.split()[0])
    if 'REFLAT=' in line:
        lat1.append(line.split()[0])
    if 'REFLONG=' in line:
        long1.append(line.split()[0])
       
# Extracts number of frequency (Nofreq)
N = nofreq1[0]
Nofreq = N[6:]
 
# Extract Latitude and Longitude 
nlat = lat1[0]         
nlong = long1[0]       
LONG = nlong[8:]
LAT = nlat[7:]

"""
  - Extracting Frequency and Real, Imaginary componets of Impedance Components
    ['FREQ',ZXXRe, ZXXImg, ZXYRe, ZXYImag,ZYXRe,ZYZImag, ZYYRe, ZYYImag] from 
    the EDI(z) file.
"""

Impedance_Components = ['ZXXR','ZXXI','ZXYR','ZXYI',
                        'ZYXR','ZYXI','ZYYR','ZYYI',
                        'FREQ']

for ii in Impedance_Components:
    for i, line in enumerate(EDIFile):    
        if '>{}'.format(ii) in line:
            impstart = i
        elif '>' in line:
            impstop = i
    Z_components1 = []
    Zs = []
         
    for j in range(impstart+1, impstop):
        raw_imp = [x for x in EDIFile[j].split('\n')]
        Z_components1.extend([l for l in raw_imp[0].split()])        
    Zs = Z_components1[:int(Nofreq)] 
    
    Zcomponents = []
    for item in Zs:
       Zcomponents.append(float(item))
       
#---------------- Write Re & Imag. Z and Freq ------------------- 
    
    with open("{}_{}.dat".format(station, ii), 'w') as text_file:    
#         text_file.write("{}\n".format(ii))   
         for j in range(len(Zs)):
             text_file.write('{:9.9s}\n'.format(str(Zcomponents[j])))
             

file_names = ['{}_FREQ.dat'.format(station),'{}_ZXXR.dat'.format(station),
           '{}_ZXXI.dat'.format(station),'{}_ZXYR.dat'.format(station),
           '{}_ZXYI.dat'.format(station),'{}_ZYXR.dat'.format(station),
           '{}_ZYXI.dat'.format(station),'{}_ZYYR.dat'.format(station),
           '{}_ZYYI.dat'.format(station)]

out = []

for files in file_names:
    filein = open(files)
    one_list = []
    one_list.append(files)
    csv_reader = csv.reader(filein, delimiter=' ')
    for row in csv_reader:
        one_list.append(row[0])
            
    out.append((n for n in one_list))# Convert the list to a generator object
    filein.close()

with open('output.txt', 'w') as fileout:
    csv_writer = csv.writer(fileout, delimiter=' ')
    for row in list(zip(*out)):
        csv_writer.writerow(row)
#  ---- Defineing variables to extract frequency and impedance components       
FREQ1 = []; ZXXR1 =[]; ZXXI1=[]; ZXYR1= []; ZXYI1=[]
ZYXR1 = []; ZYXI1 =[]; ZYYR1 = []; ZYYI1= []

for index, line in enumerate(open("output.txt")):
    if index > 0:       
        FREQ1.append(line.split()[0])
        ZXXR1.append(line.split()[1])
        ZXXI1.append(line.split()[2])
        ZXYR1.append(line.split()[3])
        ZXYI1.append(line.split()[4])
        ZYXR1.append(line.split()[5])
        ZYXI1.append(line.split()[6])
        ZYYR1.append(line.split()[7])
        ZYYI1.append(line.split()[8])
        
#  ---- converts strings to folats       
ZXXR = list(map(float, ZXXR1))
ZXXI = list(map(float, ZXXI1))
ZXYR = list(map(float, ZXYR1))
ZXYI = list(map(float, ZXYI1))
ZYXR = list(map(float, ZYXR1))
ZYXI = list(map(float, ZYXI1))
ZYYR = list(map(float, ZYYR1))
ZYYI = list(map(float, ZYYI1))
FREQ = list(map(float, FREQ1))
 

"""
   - Calculationg Complex Impedances [ZXXc,ZXYc,ZYXc,ZYYc] from the Real
     and Imaginary componets of Impedance tensor.
             ZXXc = ZXX_Re + ZXX_Imag.j
             ZXYc = ZXY_Re + ZXY_Imag.j
             ZyXc = ZYX_Re + ZYX_Imag.j
             ZXXc = ZYY_Re + ZYY_Imag.j
   - Calculates the Period (1/freq).
   - Calculates the Apparent Resisitvity and Phase using the following formulas.
     
            Phase = arctan{Z.Imag/Z.Re}, (for each component)
            Apparent Res. = 0.2 * Period*|Z,| (for each component)
            DET = 0.2*Period*|ZXXc-ZYYc - ZXYc*ZYX|, (Ranganayaki,1984; Vozoff, 
            1991)
            DET is the determinant mode (rotationally invariant)
""" 

# --- Defineing variables for complex imp. to be calculated
ZXXc = []; ZYYc = []; ZXYc = []; ZYXc = []
# ---- Defineing variables for period to be calculated
period = []
# ---- Defineing variables for resistivity to be calculated
RHOxy = []; RHOyx = []; RHOxx = []; RHOyy = []
phyYX=[]; phyXY = []; phyXX= []; phyYY= []
# ---- Defineing a variable for the derminant mode to be calculated

DET = []

for i in range(len(FREQ)):
    ZXXc.append(complex(ZXXR[i],ZXXI[i]))
    ZXYc.append(complex(ZXYR[i],ZXYI[i]))
    ZYXc.append(complex(ZYXR[i],ZYXI[i]))
    ZYYc.append(complex(ZYYR[i],ZYYI[i]))

    
    phyYX.append(math.degrees(math.atan(ZYXc[i].imag/ZYXc[i].real)))
    phyXY.append(math.degrees(math.atan(ZXYc[i].imag/ZXYc[i].real)))
    phyXX.append(math.degrees(math.atan(ZXXc[i].imag/ZXXc[i].real)))
    phyYY.append(math.degrees(math.atan(ZYYc[i].imag/ZYYc[i].real)))
    
    period.append(1/FREQ[i])
    
    RHOxy.append(0.2*period[i]*abs(ZXYc[i]**2))
    RHOyx.append(0.2*period[i]*abs(ZYXc[i]**2))
    RHOxx.append(0.2*period[i]*abs(ZXXc[i]**2))
    RHOyy.append(0.2*period[i]*abs(ZYYc[i]**2))
    
    DET.append((0.2*abs(ZXXc[i]*ZYYc[i]-ZXYc[i]*ZYXc[i]))/FREQ[i])

"""
 - Plotting Apparent Resisitvity and Phase and save the result in the working
   directory in .JPG format
"""   

plt.style.use('seaborn-white')

#------Plot Apparent resistivity xy,yx-----------------
plt.figure(figsize=(10,10))
plt.subplot(2,1,1)
plt.plot(period,RHOxy,color='blue',marker='o',linestyle =' ',label ='Res_XY')
plt.plot(period,RHOyx,color ='red',marker='o',linestyle =' ',label='Res_YX')
minval=min( min(RHOxy,RHOyx))
maxval=max(max(RHOxy,RHOyx))
plt.ylim([minval/10,maxval*10])
plt.xlim(0.5*min(period),2*max(period))
plt.xscale('log')
plt.yscale('log')
plt.xticks(size=18,weight = 'bold')
plt.yticks(size=18,weight = 'bold')
plt.legend(fontsize = 18,loc = 'best')          
plt.ylabel(s='App. Res. [$\mathbf{\Omega \cdot m}$]',size=18,weight = 'bold')
plt.grid(which ='minor',axis = 'both',alpha=.5)
plt.title('MT Station {}'.format(station),size=20,weight = 'bold')
 
#---------Plot Phase xy,yx----------------

plt.subplot(2,1,2)
plt.plot(period,phyXY,color='b',marker='s',linestyle = ' ', label = 'Phase_XY')
plt.plot(period,phyYX,color='red',marker='s',linestyle =' ',label = 'Phase_YX')
plt.xscale('log')
plt.xticks(size=15,weight = 'bold')
plt.yticks(size=15,weight = 'bold')
plt.xlabel(s='Period [s]',size=18,weight = 'bold')          
plt.ylabel(s='Phase [deg.]',size=18,weight = 'bold')
plt.grid(which ='major',axis = 'both',alpha=.5)
plt.legend(fontsize = 18,loc = 'best')
plt.savefig('{}_APPRES_PHASE.jpg'.format(station))
plt.show()

"""
#---------------Calculate Modified Impedances and Combinations----------------
#   << References>>
#   - A.Jones, mtanal.f  code from MTnet website
#   - Bahr, K. (1988). Interpretation of the magnetotelluric impedance tensor:
#     regional induction and local telluric distortion, J. Geophys. Res., 62, 
#     119–127.
#   - Chave, A. D., A. G. Jones (2012). The Magnetotelluric Method: Theory and 
#     Practice, Cambridge University Press. 
#   - Simpson, F. and Bahr, K. (2005). Practical Magnetotellurics, Cambridge  
#     University Press.
#   - Vozoff, K. (1972). The magnetotelluric method in the exploration of 
#     sedimentary  basins. Geophysics 37: 98–141.
#
# Equation (5.6) in Simpson and Bahr (2005) book defines Sum and Difference Impedances as
#                      S1 = zxx + zyy
#                      S2 = zxy + zyx
#                      D1 = zxx - zyy
#                      D2 = zxy - zyx
#
# Equation (5.20) of Simpson and Bahr (2005) show the difference between the phases of two 
# complex numbers x and y can be determined using the commutator:
#        [x,y] = xy = real(x)*aimag(y) - aimag(x)*real(y)
#
#        S1D1 = real(S1)*aimag(D1) - aimag(S1)*real(D1)
#        S2D2 = real(S2)*aimag(D2) - aimag(S2)*real(D2)
#        S1S2 = real(S1)*aimag(S2) - aimag(S1)*real(S2)
#        D1D2 = real(D1)*aimag(D2) - aimag(D1)*real(D2)
#        D1S2 = real(D1)*aimag(S2) - aimag(D1)*real(S2)
#        S1D2 = real(S1)*aimag(D2) - aimag(S1)*real(D2)
#        A  = S1D1 + S2D2
#        B  = S1S2 - D1D2
#        Cb = D1S2 - S1D2
#
"""

#  ---- Defineing variables for SUM and Difference imp. to be calculated
S1 = []; S2 = []; D1 = []; D2 = []
# ---- Defineing variable for Impedances Combinations to be calculated
S1D1 = []; S2D2 = []; S1S2 = []; D1D2 = []; D1S2=[]; S1D2 = []
# ---- Defineing variables for Impedances Combinations to be calculated
A = []; B = []; cb = []
# ---- Defineing variables for swift Skew, Phase sensitive skew and ellipticity to be calculated
Swift_Skew = []; mu = []; Ellipticity = []

for i in range(len(period)):
    S1.append(ZXXc[i]+ZYYc[i])
    S2.append(ZXYc[i]+ZYXc[i])
    D1.append(ZXXc[i]-ZYYc[i])
    D2.append(ZXYc[i]-ZYXc[i])
    S1D1.append((S1[i].real * D1[i].imag)-(S1[i].imag * D1[i].real))
    S1S2.append((S1[i].real * S2[i].imag)-(S1[i].imag * S2[i].real))    
    S2D2.append((S2[i].real * D2[i].imag)-(S2[i].imag * D2[i].real))
    D1D2.append((D1[i].real * D2[i].imag)-(D1[i].imag * D2[i].real))
    D1S2.append((D1[i].real * S2[i].imag)-(D1[i].imag * S2[i].real))
    S1D2.append((S1[i].real * D2[i].imag)-(S1[i].imag * D2[i].real))
    A.append(S1D1[i]*S2D2[i])
    B.append(S1S2[i] - D1D2[i])
    cb.append(D1S2[i] - S1D2[i])
    mu.append(((abs(cb[i]))**0.5)/abs(D2[i]))
    Swift_Skew.append(abs(S1[i])/abs(D2[i]))
    Ellipticity.append(abs(D1[i])/abs(S2[i]))
plt.figure(figsize=(10,8))

#--------Plots Swift_Skew -----------------------
plt.subplot(2,1,1)
plt.plot(period,Swift_Skew,marker = 'o',linestyle = '--',color = 'k')

plt.xscale('log')
plt.xticks(size=15,weight = 'bold')
plt.yticks(size=15,weight = 'bold')
plt.legend(fontsize = 18,loc = 'best')          
plt.ylabel(s='Swif_Skew, [kappa]',size=18,weight = 'bold')
plt.grid(which ='major',axis = 'both',alpha=.5)
plt.title('MT Station {}'.format(station),size=20,weight = 'bold')

#-----------Threshod for Phase-Sensitive Skew (mu)-------------------------- 
plt.subplot(2,1,2)
threshold_mu= []
for i in range(len(period)):
    threshold_mu.append(0.3)
    
#-----------Plots Phase-Sensitive Skew (mu)----------------
plt.plot(period,mu,marker='s',linestyle = '--',color = 'k')
#plt.plot(period,threshold_mu,color='r', linestyle = '--')

plt.xscale('log')
plt.xticks(size=15,weight = 'bold')
plt.yticks(size=15,weight = 'bold')
plt.xlabel(s='Period [s]',size=18,weight = 'bold')          
plt.ylabel(s='Bahr Skew, [mu]',size=18,weight = 'bold')
plt.grid(which ='major',axis = 'both',alpha=.5)
plt.legend(fontsize = 18,loc = 'best')
plt.savefig('{}_Kappa_mu.jpg'.format(station))
plt.show()

##-----Plots Ellipticity----------------------------------------
#plt.plot(period,Ellipticity,marker='s',linestyle = '--', 
#         label = 'e < 1--1D DIM.')
##plt.plot(period,threshold_mu,color='r', linestyle = '--')
#
#plt.xscale('log')
#plt.xticks(size=15,weight = 'bold')
#plt.yticks(size=15,weight = 'bold')
#plt.xlabel(s='Period [s]',size=18,weight = 'bold')          
#plt.ylabel(s='Ellipticity,e',size=18,weight = 'bold')
#plt.grid(which ='major',axis = 'both',alpha=.5)
#plt.legend(fontsize = 18,loc = 'best')
#
#plt.show()     


with open("{}.mt".format(station), 'w') as text_file:    
    text_file.write("# MT Station: {}\n".format(station))
    text_file.write("# No.Freq. = {}\n".format(Nofreq))   
    text_file.write("# LAT = {}, LONG = {}\n".format(LAT,LONG))
    text_file.write("\n")       
    text_file.write('   Freq[Hz]\t  RHOxy\t    RHOyx\t      phyXY\t      phyYX\t     RHOxx\t     RHOyy\t    phyXX\t    phyYY\t  Swift-Skew\tPh.-Sen.-Skew\t\n')
    text_file.write('   ----------------------------------------------------------------------------------------------------------------------------------------\n')
    for j in range(len(period)):
            text_file.write('   {:7.7s}\t {:7.7s}\t {:7.7s}\t {:7.7s}\t {:7.7s}\t {:7.7s}\t {:7.7s}\t {:7.7s}\t {:7.7s}\t {:7.7s}\t {:7.7s}\n'.format(
                         str(FREQ[j]),str(RHOxy[j]),str(RHOyx[j]),str(phyXY[j]),str(phyYX[j]),str(RHOxx[j]),str(RHOyy[j]),str(phyXX[j]),str(phyYY[j]),str(Swift_Skew[j]),str(mu[j])))

#To remove created .dat files in the working directory

for ii in Impedance_Components:
    os.remove("{}\{}_{}.dat".format(name_of_folder,station,ii))
os.remove("output.txt")
