
# This Python (v2.7) Program cacluates Swift-Skew and Phase Sen. Skew from all edi files saved in a folder. 

# Author: Biruk A. Cherkose, @2019, Email: 201990083@uaeu.ac.ae

#-----------------------------INPUT from the USER------------------------------
# write your folder location where the input EDI files are saved
# EDI files should be in Impedance format

name_of_folder = r'.'

#-------------------------------------------------------------------------------


import math
import os
import glob
import csv
import matplotlib.pyplot as plt



input_files = glob.glob(os.path.join(name_of_folder, '*.edi'))

for indexs in input_files:
    station = '{}'.format(indexs)
    f = open("{}".format(station))
    EDIFile = f.readlines()

    nofreq1 = []
    lat1=[]
    long1 = []

    for index, line in enumerate(open("{}".format(station))):
        if 'NFREQ=' in line:
            nofreq1.append(line.split()[0])
        if 'REFLAT=' in line:
            lat1.append(line.split()[0])
        if 'REFLONG=' in line:
            long1.append(line.split()[0])

    N = nofreq1[0]
    Nofreq = N[6:]
    nlat = lat1[0]
    LAT = nlat[7:]
    nlong = long1[0]
    LONG = nlong[8:]
    """
     -------Extracting Frequencies from the EDI files-------------------------
    """
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
    """
       Extracting Real and Imaginary componets of Impedance Components
       [ZXXRe, ZXXImg, ZXYRe, ZXYImag,ZYXRe,ZYZImag, ZYYRe, ZYYImag]
    """        
    for i, line in enumerate(EDIFile):
        if '>ZXXR' in line:
            impstart = i
        elif '>' in line:
            impstop = i
    ZXXR1 = []
    for j in range(impstart+1, impstop):
        raw_imp = [k for k in EDIFile[j].split('\n')]
        ZXXR1.extend([l for l in raw_imp[0].split()])

    ZXXReal = ZXXR1[:int(Nofreq)]
    for i, line in enumerate(EDIFile):
        if '>ZXXI' in line:
            impstart = i
        elif '>' in line:
            impstop = i
    ZXXI1 = []
    for j in range(impstart+1, impstop):
        raw_imp = [k for k in EDIFile[j].split('\n')]
        ZXXI1.extend([l for l in raw_imp[0].split()])

    ZXXImag = ZXXI1[:int(Nofreq)]

    for i, line in enumerate(EDIFile):
        if '>ZXYR' in line:
            impstart = i
        elif '>' in line:
            impstop = i
    ZXYR1 = []
    for j in range(impstart+1, impstop):
        raw_imp = [k for k in EDIFile[j].split('\n')]
        ZXYR1.extend([l for l in raw_imp[0].split()])

    ZXYReal = ZXYR1[:int(Nofreq)]
    for i, line in enumerate(EDIFile):
        if '>ZXYI' in line:
            impstart = i
        elif '>' in line:
            impstop = i
    ZXYI1 = []

    for j in range(impstart+1, impstop):
        raw_imp = [k for k in EDIFile[j].split('\n')]
        ZXYI1.extend([l for l in raw_imp[0].split()])

    ZXYImag = ZXYI1[:int(Nofreq)]
    for i, line in enumerate(EDIFile):
        if '>ZYXR' in line:
            impstart = i
        elif '>' in line:
            impstop = i
    ZYXR1 = []
    for j in range(impstart+1, impstop):
        raw_imp = [k for k in EDIFile[j].split('\n')]
        ZYXR1.extend([l for l in raw_imp[0].split()])
    ZYXReal = ZYXR1[:int(Nofreq)]
    for i, line in enumerate(EDIFile):
        if '>ZYXI' in line:
            impstart = i
        elif '>' in line:
            impstop = i
    ZYXI1 = []
    for j in range(impstart+1, impstop):
        raw_imp = [k for k in EDIFile[j].split('\n')]
        ZYXI1.extend([l for l in raw_imp[0].split()])
    ZYXImag = ZYXI1[:int(Nofreq)]
    for i, line in enumerate(EDIFile):
        if '>ZYYR' in line:
            impstart = i
        elif '>' in line:
            impstop = i
    ZYYR1 = []
    for j in range(impstart+1, impstop):
        raw_imp = [k for k in EDIFile[j].split('\n')]
        ZYYR1.extend([l for l in raw_imp[0].split()])
    ZYYReal = ZYYR1[:int(Nofreq)]
    for i, line in enumerate(EDIFile):
        if '>ZYYI' in line:
            impstart = i
        elif '>' in line:
            impstop = i
    ZYYI1 = []
    for j in range(impstart+1, impstop):
        raw_imp = [k for k in EDIFile[j].split('\n')]
        ZYYI1.extend([l for l in raw_imp[0].split()])
    ZYYImag = ZYYI1[:int(Nofreq)]
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
    """
     - Calculating Complex Impedances [ZXXc,ZXYc,ZYXc,ZYYc] from the Real
       and Imaginary componets of Impedance tensor.
             ZXXc = ZXX_Re + ZXX_Imag.j
             ZXYc = ZXY_Re + ZXY_Imag.j
             ZyXc = ZYX_Re + ZYX_Imag.j
             ZXXc = ZYY_Re + ZYY_Imag.j
     - Calculates the Period (1/freq).
     - Calculates the Apparent Resisitvity and Phase from each  EDI files in the folder.

             Phase = arctan{Z.Imag/Z.Re}, (for each component)
             Apparent Res. = 0.2 * Period*|Z| (for each component)
             DET = 0.2*Period*|ZXXc-ZYYc - ZXYc*ZYX|, (Ranganayaki,1984; Vozoff, 1991)

    """
    zxx = []; zyy = []; zxy = []; zyx = []
    ZXXc = []; ZYYc = []; ZXYc = []; ZYXc = []

    period = []

    RHOxy = []; RHOyx = []; RHOxx = []; RHOyy = []
    phyYX=[]; phyXY = []; phyXX= []; phyYY= []
            #  ---- Defineing variables for SUM and Difference to be calculated
    S1 = []; S2 = []; D1 = []; D2 = []
        # ---- Defineing variable for ImpedancesCombinations to be calculated
    S1D1 = []; S2D2 = []; S1S2 = []; D1D2 = []; D1S2=[]; S1D2 = []
        # ---- Defineing variables for Impedances Combinations to be calculated
    A = []; B = []; cb = []
        # ---- Defineing variables for swift Skew, Phase sensitive skew and ellipticity to be calculated
    Swift_Skew = []; mu = []; Ellipticity = [] 


    DET = []
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

        DET.append((0.2*abs(ZXXc[i]*ZYYc[i]-ZXYc[i]*ZYXc[i]))/freq[i])
    for i in range(len(period)):


    # -------------------Calculating Modified Impedances and Combinations-----------

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

    #---------Calculating swift Skew, Phase sensitive skew (mu) and ellipticity----
        mu.append(((abs(cb[i]))**0.5)/abs(D2[i]))
        Swift_Skew.append(abs(S1[i])/abs(D2[i]))
        Ellipticity.append(abs(D1[i])/abs(S2[i]))        

    with open("{}_mu_Kappa_Out.dat".format(station), 'w') as text_file:
        for j in range(len(freq)):
                text_file.write('{:7.7s}\t {:7.7s}\t {:7.7s}\n'.format(str(period[j]),
                             str(Swift_Skew[j]),str(mu[j])))
                
                
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
output_filename = '{}\Combine_Dim.txt'.format(name_of_folder)

#-------------------------------------------------------------------------------
input_files = glob.glob(os.path.join(name_of_folder, '*.dat'))
with open(os.path.join(name_of_folder, output_filename), 'w') as file_out:
     headers_read = False
     for input_file in input_files:
         if input_file == os.path.join(name_of_folder, output_filename):
# If the output file is in the list of input files, ignore it.
             continue
         with open(input_file, 'r') as fin:
              reader = csv.reader(fin)
              if not headers_read:
# Read column headers just once
                  headers = reader.next()[0].split()
                  headers = headers[0:2] + [headers[2]]
                  file_out.write("\t".join(headers + ['\n'])) # Zero based indexing.
                  headers_read = True

#              else:
#                   _ = reader.next() # Ignore header row.
              for line in reader:
                  if line: # Ignore blank lines.
                     line_out = line[0].split()
                     file_out.write("\t".join(line_out[0:2] + [line_out[2]] + ['\n']))                
                
fileopen = r"Combine_Dim.txt"

pp1 = []
ss1 = [] 
mu11 = []  
              
for index, line in enumerate(open(fileopen)):
    if index >=0:
        pp1.append(line.split()[0])
        ss1.append(line.split()[1])
        mu11.append(line.split()[2])
#print mu11

period_new = []
swift_Skew_new = []
Phase_Sen_Skew_new = []

for item in pp1:                  
    period_new.append(float(item))          
                        
for item in ss1:                 
    swift_Skew_new.append(float(item))
for item in mu11:                  
    Phase_Sen_Skew_new.append(float(item))  
#print period_new
thss = []
for num in range(len(period_new)):
    thss.append(0.2)
thpss = []
for num in range(len(period_new)):
    thpss.append(0.3)
    
    
# _plot Result 
             
plt.style.use('seaborn-white')
plt.figure(figsize=(10,12))
plt.subplot(211)
plt.plot(period_new,swift_Skew_new,
         label='Swift_Skew',
         marker = 'o',
         linestyle = ' ')
plt.plot(period_new,thss,
         color = 'r',
         linestyle = '-',
         lw = 0.3)
#plt.box()
plt.xscale('log')
plt.xticks(size=18)
plt.yticks(size=18)
plt.legend(fontsize=17)
#minval=min( min(obs,pre))
#maxval=max(max(obs,pre))
plt.ylim(0,1.5)
#plt.xlim(2*max(freq),0.5*min(freq))
plt.grid()
#plt.title(Title2, size=20,weight = 'bold')
#plt.xlabel(s='Period [s]',
#           size=18,
#           weight = 'bold')
plt.ylabel(s='Swift_Skew (Kappa)',
           size=18,
           weight = 'bold')



#++++++++++++++++++++++++++++++++++++++
plt.subplot(212)
plt.plot(period_new,Phase_Sen_Skew_new,
         label='mu',
         color = 'g',
         marker = 'o',
         linestyle = ' ')
plt.plot(period_new,thpss,
         color = 'r',
         linestyle = '-',
         lw = 0.3)

#plt.box()
plt.xscale('log')
plt.xticks(size=18)
plt.yticks(size=18)
plt.legend(fontsize=17)
#minval=min( min(obs,pre))
#maxval=max(max(obs,pre))
#plt.ylim([minval/10,maxval*10])
#plt.xlim(2*max(freq),0.5*min(freq))
plt.grid()
#plt.title(Title2, size=20,weight = 'bold')
plt.xlabel(s='Period [s]',
           size=18,
           weight = 'bold')
plt.ylabel(s="Phase Sensitive Skew (mu)",
           size=18,
           weight = 'bold')
plt.savefig("All_mu_Kappa.jpeg")

plt.show()
