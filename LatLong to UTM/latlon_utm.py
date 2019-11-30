#==============================================================================
#              ***  latlon_utm.py   ***
#
#    This Python (V2.7) program converts lat&long (degree, min., sec.) in the edi files
#    to utm(meters, WGS84) and produce one file ("Coordinate_Lat-Long_UTM.dat") that contains 
#    utm(X&Y) from all edi files in a directory. 
#    
#    Search LAT and LONG in  edi files in a folder and convert them into utm
#    (WGS84). 
#    
#    HOW TO USE THE PROGRAM
#    INPUT from the USER
#        * The user input the Project Area
#        * Path to edi files
#    OUTPUT
#          * A file that contains the lat/long and utm(Easting & Northing) from all EDI files
#          * File produced in the working directory(Coordinate_Lat-Long_UTM.dat)
#
#  <<Refrences>>
#      - MTpy: A python tool for MT data analysis and processing
#      - A fuction from a script, conversions.py (from MTpy) is used, @UofA, 2013 (LK). 
#        for coversion purpose lat/long to utm (WGS84).
#               
#      
#
# Prepared by: Biruk Abera (June, 2019)
#==============================================================================

import os
import glob
import conversions as utm2ll

#----------------------USER Inputs Section --------------------
Project_Area = "Project Name" # Write Project name
name_of_folder = r'.' # Write path of EDI folder
#---------------------------------------------------------------



input_files = glob.glob(os.path.join(name_of_folder, '*.edi'))

for indexs in input_files:

    station = '{}'.format(indexs)
    
    f = open("{}".format(station))
    EDIFile = f.readlines()

    lat1=[]
    long1 = []
    station_point1 = []
    for index, line in enumerate(open("{}".format(station))):
        if 'REFLAT=' in line:
            lat1.append(line.split()[0])
        if 'REFLONG=' in line:
            long1.append(line.split()[0])
        if 'DATAID=' in line:
            station_point1.append(line.split()[0])
            
    sp = station_point1[0]
    station_point = sp[7:]
    nlat = lat1[0]
    LAT = nlat[7:]
    nlong = long1[0]
    LONG = nlong[8:]

    with open("coordinate.txt", "a+") as textfile:
         textfile.write("{}\t{}\t{}\n".format(LAT, LONG,station_point))

ff = open(r"coordinate.txt", "r")     
la1=[]
lo1 = []
spedi = []

for index, line in enumerate(ff):
    la1.append(line.split()[0])
    lo1.append(line.split()[1])
    spedi.append(line.split()[2])

ff.close()
#----------------LAT------------
deg1 = []; m1 = []; s1 = []; xx = []
for iv in la1:
    deg1.append(float(iv[:iv.index(":")]))
    m1.append(float(iv[iv.index(":")+1:-iv[::-1].index(":")-1])/60)
    s1.append(float(iv[-iv[::-1].index(":"):-1])/3600)
#    xx.append(float(deg1[iv]) + float(m1[iv]) + float(s1[iv]))
for iii in range(len(la1)):
    xx.append(deg1[iii]+ m1[iii] + s1[iii])
#-----------LONG    
deg12 = []; m12 = []; s12 = []; yy = []
for ix in lo1:
    deg12.append(float(ix[:ix.index(":")]))
    m12.append(float(ix[ix.index(":")+1:-ix[::-1].index(":")-1])/60)
    s12.append(float(ix[-ix[::-1].index(":"):-1])/3600)
#    xx.append(float(deg1[iv]) + float(m1[iv]) + float(s1[iv]))
for ixx in range(len(lo1)):
    yy.append(deg12[ixx]+ m12[ixx] + s12[ixx])

for jj in range(len(input_files)):    
    UTM_EN= utm2ll.LLtoUTM(23, xx[jj], yy[jj])

    
    with open("utm.txt", "a+") as textfile2:
        textfile2.write("{}\n".format(list(UTM_EN)))
        
with open('utm.txt', 'r') as file3:
  filedata = file3.read()

# Replace the target string
filedata = filedata.replace(']', '0')
filedata = filedata.replace('[', '0')
filedata = filedata.replace(',', '0')
with open('UTMcon.txt', 'w') as file3:
  file3.write(filedata)
    
         
UTM_northing=[]
UTM_easting = [] 
fff = open("UTMcon.txt", "r")
        
for index, line in enumerate(fff):
    UTM_northing.append(line.split()[2])
    UTM_easting.append(line.split()[1])
UTM_E = []; UTM_N = []
for items in UTM_easting:
    UTM_E.append(float(items))
for items in UTM_northing:
    UTM_N.append(float(items))

with open("Coordinate_Lat-Long_UTM.dat", "w") as textfile000:
    textfile000.write(" # Project Area: {}\n".format(Project_Area))
    textfile000.write(" # Produced by latlon_utm.py(@#Biruk)\n".format(Project_Area))
    textfile000.write('\n')    
    textfile000.write(' Station      LAT      LONG       EASTING(WGS84)   NORTHING(WGS84) \n')      
    textfile000.write(' ---------------------------------------------------------------\n')
    textfile000.write('\n')    
    for ss in range(len(input_files)):
        textfile000.write(" {}\t{:9.9s}\t{:9.9s}\t{:9.9s}\t{:9.9s}\t\n".format(str(spedi[ss]),
                          str(xx[ss]),str(yy[ss]),str(UTM_E[ss]),str(UTM_N[ss])))       
fff.close()
    
os.remove("coordinate.txt")  
os.remove("utm.txt")    
os.remove("UTMcon.txt")  

#for file in os.listdir('.\EDI'):
#    print file