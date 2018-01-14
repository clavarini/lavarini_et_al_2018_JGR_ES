##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## Pebble abrasion model 
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## This model calculates the contribution of each lithological unity to the  
## proportion of heavy minerals present in modern fluvial sands.
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
## We do this by calculating the following parameters from every lithological 
## unity: exposure area, specific heavy mineral concentration in the original 
## rock, initial bedload to suspended load ratio, rock erodibility to abrasion, 
## erosion rate and sediment travel distance.
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## The model reads four database files originally from ArcGIS where exposure 
## area, lithological unities, elevation, flow length and flow accumulation 
## with their spatial reference (X, Y) are added in the model. 
## The remaining parameters (erodibility, load ratio, mineral concentration, 
## and erosion rate) are user-dependent. 
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## As an example, this model is tested to figure out the influence of
## abrasion on the zircon population of sands in the Marsyandi watershed, Nepal.
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## Authors:
## Mikael Attal, Chrystiann Lavarini & Carlos A. da Costa Filho
## The University of Edinburgh, UK
## January of 2017
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# Importing modules to perform calculations.

import numpy as np

###############################################################################
#                                                                             #
#                             READ IN THE DATA                                #
#                                                                             #
###############################################################################

#   Object attributes are computed from file directory and DBF filenames      #
#                    containing the model dataset.                            #


DataDirectory = r'D:/Chrystiann Lavarini/Universidade/Doutorado/PhD Project/MODEL/' 

# These files were extracted through ModelBuilder in ArcGIS.                  #

FileName = 'MARS_90_DEM_PT.dbf' # They represent elevation.
FileName1 = 'F_ACC_MARS_90_PT.dbf' # They represent flow accumulation.
FileName2 = 'EXT_GEO_MARS_K_BASIN.dbf'  # They represent lithological unity with extended Formation II/III.
#FileName2 = 'GEO_MARS_ORI_90_.dbf'  # They represent lithological unity with original Formation II/III.
FileName3 = 'F_LEN_MARS_90_PT.dbf' # They represent flow length.

# Creates series of arrays to perform calculations.                           #
                                                                                        
PC_zir_ALL = [] 
bedloadratio_ALL = [] 
PC_sus_ALL = [] 
PC_bed_ALL = [] 
PC_pixels_ALL = [] 

###############################################################################
#                                                                             #                     
#                             Using dbfread                                   #
#                                                                             #                     
###############################################################################

from dbfread import DBF

# Spatial coordinantes (x, y) and elevation (z) are collected (from FileName).# 

dbffile = DBF(DataDirectory + FileName)
x = []
y = []
z = []

for record in dbffile:
    x.append(record['x'])
    y.append(record['y'])
    z.append(record['grid_code'])

#   Area (area) of every lithological unity is collected (from FileName1).    #

dbffile = DBF(DataDirectory + FileName1)
area = []
for record in dbffile:
    area.append(record['grid_code'])

#        Every lithological unity is collected (from FileName2).              #

dbffile = DBF(DataDirectory + FileName2)
litho = []
for record in dbffile:
    litho.append(record['RASTERVALU'])
    

LITHO = np.array(litho)    
X = np.array(x)
Y = np.array(y)

###############################################################################        
for n,i in enumerate(LITHO): # A subsitution is made here because some        #
    if i== -9999:            # discrepancies (-9999 error) occured when ArcGIS#
      LITHO[n]=1             # mapped the lithological unities. The error     #
                             # should be fixed by lithology 1 (Tethys unity). #
###############################################################################        
      
# Flow length is collected (from FileName3)                                   #

dbffile = DBF(DataDirectory + FileName3)
length = []
for record in dbffile:
    length.append(record['grid_code'])

###############################################################################        
#                                                                             #
#                  Input user-dependent parameters                            #
#                                                                             #
###############################################################################        
   
#Set parameters for model:
    
pixel_size=90 # Pixel size in meters, used to calculate volumes but not     #
                # essential for percentages provided at the end.              #
eros=0.001 # uniform erosion rate, in m/yr, used to calculate volumes but not #
           # essential for percentages provided at the end.                   #

###############################################################################        
#  future development: spatially variable erosion based on published data     #
###############################################################################        

n_lithos = len(set(LITHO))
litho_names = []

for i in LITHO:
    if LITHO[i] ==1:
        litho_names.append('Tethyan series')
    elif LITHO[i] ==2:
        litho_names.append('Formation II/III')
    elif LITHO[i] ==3:
        litho_names.append('Manaslu granite')
    elif LITHO[i] ==4:
        litho_names.append('Formation I')
    else:
        litho_names.append('Lesser Hyamalaia')

###############################################################################
# Abrasion coefficient from Attal & Lave (2006)   for the ORIGINAL 
###############################################################################
   
#erod1=4.30  #      erodibility coefficient of litho 1, in % mass loss/km    #
#erod2=0.40  #           erodibility coefficient of litho 2.                 #
#erod3=0.40  #             erodibility coefficient of litho 3.               #
#erod4=1.24  #             erodibility coefficient of litho 4.               #
#erod5=9.42 #            erodibility coefficient of litho 5.                #

###############################################################################
# Abrasion coefficient from Attal & Lave (2006)   for the EXTENDED
###############################################################################

#erod1=4.3  #      erodibility coefficient of litho 1, in % mass loss/km    #
#erod2=0.4  #      erodibility coefficient of litho 1, in % mass loss/km    #
#erod3=0.4 #             erodibility coefficient of litho 3.               #
#erod4=1.23 #             erodibility coefficient of litho 4.               #
#erod5=9.42  #            erodibility coefficient of litho 5.                #
   
# This averages were calculated from a weighted arithmetic mean of Attal and Lave (2006).

erod1=0  #      erodibility coefficient of litho 1, in % mass loss/km    #
erod2=0  #           erodibility coefficient of litho 2.                 #
erod3=0 #             erodibility coefficient of litho 3.               #
erod4=0  #             erodibility coefficient of litho 4.               #
erod5=0 #            erodibility coefficient of litho 5.                #

erod=[erod1, erod2, erod3, erod4, erod5]

fert1=0.1        #    fraction of rock type 1 that is zircon ('fertility') #
fert2=0.1       #           fraction of rock type 2 that is zircon        #
fert3=0      #           fraction of rock type 3 that is zircon        #
fert4=0.1       #           fraction of rock type 4 that is zircon        #
fert5=0.1       #           fraction of rock type 5 that is zircon        #

###############################################################################
#  Data of Amidon et al. (2005b) where grain/gram to GEOL_1  #
###############################################################################

#fert1=5.6        #    fraction of rock type 1 that is zircon ('fertility') #
#fert2=0.83        #           fraction of rock type 2 that is zircon        #
#fert3=0.00       #           fraction of rock type 3 that is zircon        #
#fert4=6.48       #           fraction of rock type 4 that is zircon        #
#fert5=2.13       #           fraction of rock type 5 that is zircon        #

###############################################################################
#  Data of Amidon et al. (2005b) where grain/gram to LITHO_POINT  #
###############################################################################

#fert1=5.6        #    fraction of rock type 1 that is zircon ('fertility') #
#fert2=3.1        #           fraction of rock type 2 that is zircon        #
#fert3=0.00       #           fraction of rock type 3 that is zircon        #
#fert4=6.48       #           fraction of rock type 4 that is zircon        #
#fert5=2.13       #           fraction of rock type 5 that is zircon        #

###############################################################################
#  Data of Amidon et al. (2005b) where data was converted to weight percent   #
###############################################################################

#fert1=0.36        #    fraction of rock type 1 that is zircon ('fertility') #
#fert2=0.07        #           fraction of rock type 2 that is zircon        #
#fert3=0.00       #           fraction of rock type 3 that is zircon        #
#fert4=0.465       #           fraction of rock type 4 that is zircon        #
#fert5=0.145       #           fraction of rock type 5 that is zircon        #

###############################################################################
#  Data of Amidon et al. (2005b) where data was kept as ratios  for litho1    #
###############################################################################

#fert1=0.57        #    fraction of rock type 1 that is zircon ('fertility') #
#fert2=0.31        #           fraction of rock type 2 that is zircon        #
#fert3=0.00000       #           fraction of rock type 3 that is zircon        #
#fert4=0.81       #           fraction of rock type 4 that is zircon        #
#fert5=0.32       #           fraction of rock type 5 that is zircon        #

###############################################################################
#  Data of Amidon et al. (2005b) where data was kept as ratios  for geol_1    #
###############################################################################

#fert1=0.57        #    fraction of rock type 1 that is zircon ('fertility') #
#fert2=0.1        #           fraction of rock type 2 that is zircon        #
#fert3=0.00000       #           fraction of rock type 3 that is zircon        #
#fert4=0.81       #           fraction of rock type 4 that is zircon        #
#fert5=0.32       #           fraction of rock type 5 that is zircon        #
###############################################################################
#       Unplished data of France-Lanord (2005b)on weight percent              #
###############################################################################

#fert1=0.54        #    fraction of rock type 1 that is zircon ('fertility') #
#fert2=0.22       #           fraction of rock type 2 that is zircon        #
#fert3=0.00       #           fraction of rock type 3 that is zircon        #
#fert4=0.57       #           fraction of rock type 4 that is zircon        #
#fert5=0.165       #           fraction of rock type 5 that is zircon        #


fert=[fert1, fert2, fert3, fert4, fert5]


source_1mm1= 0.25 # fraction of source material finer than 1 mm for rock type 1
source_1mm2= 0.25 # fraction of source material finer than 1 mm for rock type 2
source_1mm3= 0.25 # fraction of source material finer than 1 mm for rock type 3
source_1mm4= 0.25 # fraction of source material finer than 1 mm for rock type 4
source_1mm5= 0.25 # fraction of source material finer than 1 mm for rock type 5

source_1mm=[source_1mm1, source_1mm2, source_1mm3, source_1mm4, source_1mm5]
    # future development: spatially variable, e.g., glacial vs landslide.
    
###############################################################################        
#                                                                             #
#               Iteration of files extracted from ArcGIS                      #
#                                                                             #
###############################################################################          

max_area=0.0; outletID=0;
n_lines = len(x)   # get the number of pixels (=number of data)               #
for i in range(0, n_lines):
    if area[i] > max_area:  #  finding the outlet node and basin total area   #
        max_area=area[i]; outletID=i 
        
###############################################################################          
#               Printing and iterating the model outcomes                     #
###############################################################################
        
print() 
print( "Marsyandi outlet") # Write a code that extract the same names that Amidon (2005b) gave to his samples. 
print()            
print( "Outlet node ID is "+str(outletID)+" at x="+str(x[outletID])+" m and \
y="+str(y[outletID])+" m; drainage area is "+str(area[outletID])+" pixels.")

#                Setting the length of all nodes to outlet node               #

length2outlet=[]
for i in range(len(length)):
    length2outlet.append(length[i]-length[outletID])
    if length[i]-length[outletID]<0.0:
      print( "Error: all nodes should be upstream of outlet") 
                                                                      
#                   Make the calculations for each lithology                  #
#                           initialise variables                              #
#             Volume arrays that will contain data for each lithology         #
#                total volumes combining all lithos ALL PER YEAR              #
                                                                      
npixels=[]; vol_bed=[]; vol_sus=[]; vol_bedsus=[]; vol_zir=[]; 
vol_bedT=0; vol_susT=0; vol_bedsusT=0; vol_zirT=0; 

for i in range(1, n_lithos+1):    # proceed for each litho
    count=0; bedload=0; suspended=0; 
    print( 'Proceeding with litho number '+str(i))
    for j in range(len(x)):     
        if litho[j]==i:         # find the nodes that are of litho type i
            # add contribution of each pixel to bedload:            
            bedload=bedload + (1-source_1mm[i-1])*eros*pixel_size*pixel_size*\
                np.exp(-length2outlet[j]/1000*erod[i-1]/100)
            # add contribution of each pixel to suspended load:
            suspended=suspended + source_1mm[i-1]*eros*pixel_size*pixel_size +\
                (1-source_1mm[i-1])*eros*pixel_size*pixel_size*\
                (1-np.exp(-length2outlet[j]/1000*erod[i-1]/100))
            # count the number of pixels of lithology i:
            count = count + 1
    
#     Once gone through the list of pixels, record the values in arrays:      #
#     First value in each aray is for litho 1, second for litho 2, etc..      #
#         (though in python it's 0th value, 1st value, etc)                   #
        
    npixels.append(count)   # arrays contains number of pixels for each litho #
    vol_bed.append(bedload) # volume of bedload for each litho (m3/yr)        #
    vol_sus.append(suspended) # volume of suspended load for litho (m3/yr)    #
    total_volumeforlitho=bedload+suspended # contains total volume of sediment#   
    vol_bedsus.append(total_volumeforlitho)  # for each litho (m3/yr)         #
    total_zirconforlitho=suspended*fert[i-1] # contains volume of zircon in   #
    vol_zir.append(total_zirconforlitho)     # sand for each litho (m3/yr)    # 
   
#                           Check for mass balance                            #
    check=suspended+bedload-(count*eros*pixel_size*pixel_size)
    if abs(check)>1:
        print('Mass balance problem, ' +str(check)+ ' m3 lost')

#            Computing total volumes out of catchment (at outlet)             #

vol_bedT=sum(vol_bed)
vol_susT=sum(vol_sus)
vol_bedsusT=sum(vol_bedsus)
vol_zirT=sum(vol_zir)

#                           Check for mass balance                            #

check=len(x)-sum(npixels)
if abs(check)>0:
    print('Calculation problem, ' +str(check)+ ' pixels lost')
check=vol_bedsusT-len(x)*eros*pixel_size*pixel_size
if abs(check)>0:
    print('Mass balance problem, ' +str(check)+ ' m3 lost')

#                            Calculate bedload ratio                          #

bedloadratio = vol_bedT/vol_bedsusT*100
print('At the outlet considered, bedload represents ' \
+ str("%.1f" % bedloadratio) + ' % of the sediment.')

#              PERCENT arrays that contain data for each litho                #

PC_pixels=[]; PC_bed=[]; PC_sus=[]; PC_bedsus=[]; PC_zir=[]; 
for i in range(len(npixels)):
    PC_pixels.append(100*float(npixels[i])/float(len(x)))
    PC_bed.append(100*vol_bed[i]/vol_bedT)
    PC_sus.append(100*vol_sus[i]/vol_susT)
    PC_bedsus.append(100*vol_bedsus[i]/vol_bedsusT)
    PC_zir.append(100*vol_zir[i]/vol_zirT)    
        
#                           Displaying the results                            #
        
for i in range(len(npixels)):
    print('LITHOLOGY '+str(i+1)+':')
    print('- makes up '+str("%.1f" % PC_pixels[i])+'% of catchment.')
    print('- makes up '+str("%.1f" % PC_bed[i])+'% of bedload, '\
    +str("%.1f" % PC_sus[i])+'% of suspended load and '\
    +str("%.1f" % PC_bedsus[i])+'% of total load.')
    print( '- makes up '+str("%.1f" % PC_zir[i])+'% of zircons found in sands.')
    



###############################################################################
#                                                                             #
#                           THE END OF THE MODEL                              #
#                                                                             #
###############################################################################

def cost_function_mistmatch(erod):
    eros  = [1, 1, 1, 1, 1]
                                                                                      
    PC_zir_ALL = [] 
    bedloadratio_ALL = [] 
    PC_sus_ALL = [] 
    PC_bed_ALL = [] 
    PC_pixels_ALL = [] 
        
#    fert1=0.57        #    fraction of rock type 1 that is zircon ('fertility') #
#    fert2=0.31        #           fraction of rock type 2 that is zircon        #
#    fert3=0.00000       #           fraction of rock type 3 that is zircon        #
#    fert4=0.81       #           fraction of rock type 4 that is zircon        #
#    fert5=0.32       #           fraction of rock type 5 that is zircon        #    

    fert1=0.1        #    fraction of rock type 1 that is zircon ('fertility') #
    fert2=0.1       #           fraction of rock type 2 that is zircon        #
    fert3=0      #           fraction of rock type 3 that is zircon        #
    fert4=0.1       #           fraction of rock type 4 that is zircon        #
    fert5=0.1       #           fraction of rock type 5 that is zircon        #
#    
    fert=[fert1, fert2, fert3, fert4, fert5]
    
#    erod1=4.30  #      erodibility coefficient of litho 1, in % mass loss/km    #
#    erod2=3.06  #      erodibility coefficient of litho 1, in % mass loss/km    #
#    erod3=0.40  #             erodibility coefficient of litho 3.               #
#    erod4=1.40  #             erodibility coefficient of litho 4.               #
#    erod5=9.42 #            erodibility coefficient of litho 5.                #

#    erod1=0  #      erodibility coefficient of litho 1, in % mass loss/km    #
#    erod2=0  #      erodibility coefficient of litho 1, in % mass loss/km    #
#    erod3=0  #             erodibility coefficient of litho 3.               #
#    erod4=0  #             erodibility coefficient of litho 4.               #
#    erod5=0 #            erodibility coefficient of litho 5.                #
#
#    erod=[erod1, erod2, erod3, erod4, erod5]

    source_1mm1=0.25 # fraction of source material finer than 1 mm for rock type 1
    source_1mm2=0.25 # fraction of source material finer than 1 mm for rock type 2
    source_1mm3=0.25 # fraction of source material finer than 1 mm for rock type 3
    source_1mm4=0.25 # fraction of source material finer than 1 mm for rock type 4
    source_1mm5=0.25 # fraction of source material finer than 1 mm for rock type 5
    
    source_1mm=[source_1mm1, source_1mm2, source_1mm3, source_1mm4, source_1mm5]
        # future development: spatially variable, e.g., glacial vs landslide.
        
    ###############################################################################        
    #                                                                             #
    #               Iteration of files extracted from ArcGIS                      #
    #                                                                             #
    ###############################################################################          
    
    max_area=0.0; outletID=0;
    n_lines = len(x)   # get the number of pixels (=number of data)               #
    for i in range(0, n_lines):
        if area[i] > max_area:  #  finding the outlet node and basin total area   #
            max_area=area[i]; outletID=i 
            
    ###############################################################################          
    #               Printing and iterating the model outcomes                     #
    ###############################################################################
            
#    print() 
#    print( "Marsyandi outlet") # Write a code that extract the same names that Amidon (2005b) gave to his samples. 
#    print()            
#    print( "Outlet node ID is "+str(outletID)+" at x="+str(x[outletID])+" m and \
#    y="+str(y[outletID])+" m; drainage area is "+str(area[outletID])+" pixels.")
    
    #                Setting the length of all nodes to outlet node               #
    
    length2outlet=[]
    for i in range(len(length)):
        length2outlet.append(length[i]-length[outletID])
        if length[i]-length[outletID]<0.0:
          print( "Error: all nodes should be upstream of outlet") 
                                                                          
    #                   Make the calculations for each lithology                  #
    #                           initialise variables                              #
    #             Volume arrays that will contain data for each lithology         #
    #                total volumes combining all lithos ALL PER YEAR              #
                                                                          
    npixels=[]; vol_bed=[]; vol_sus=[]; vol_bedsus=[]; vol_zir=[]; 
    vol_bedT=0; vol_susT=0; vol_bedsusT=0; vol_zirT=0; 
    
    for i in range(1, n_lithos+1):    # proceed for each litho
        count=0; bedload=0; suspended=0; 
#        print( 'Proceeding with litho number '+str(i))
        for j in range(len(x)):     
            if litho[j]==i:         # find the nodes that are of litho type i
                # add contribution of each pixel to bedload:            
                bedload=bedload + (1-source_1mm[i-1])*eros[i-1]*pixel_size*pixel_size*\
                    np.exp(-length2outlet[j]/1000*erod[i-1]/100)
                # add contribution of each pixel to suspended load:
                suspended=suspended + source_1mm[i-1]*eros[i-1]*pixel_size*pixel_size +\
                    (1-source_1mm[i-1])*eros[i-1]*pixel_size*pixel_size*\
                    (1-np.exp(-length2outlet[j]/1000*erod[i-1]/100))
                # count the number of pixels of lithology i:
                count = count + 1
        
    #     Once gone through the list of pixels, record the values in arrays:      #
    #     First value in each aray is for litho 1, second for litho 2, etc..      #
    #         (though in python it's 0th value, 1st value, etc)                   #
            
        npixels.append(count)   # arrays contains number of pixels for each litho #
        vol_bed.append(bedload) # volume of bedload for each litho (m3/yr)        #
        vol_sus.append(suspended) # volume of suspended load for litho (m3/yr)    #
        total_volumeforlitho=bedload+suspended # contains total volume of sediment#   
        vol_bedsus.append(total_volumeforlitho)  # for each litho (m3/yr)         #
        total_zirconforlitho=suspended*fert[i-1] # contains volume of zircon in   #
        vol_zir.append(total_zirconforlitho)     # sand for each litho (m3/yr)    # 
       
    #                           Check for mass balance                            #
        check=suspended+bedload-(count*eros[i-1]*pixel_size*pixel_size)
#        if abs(check)>1:
#            print('Mass balance problem, ' +str(check)+ ' m3 lost')
#    
    #            Computing total volumes out of catchment (at outlet)             #
    
    vol_bedT=sum(vol_bed)
    vol_susT=sum(vol_sus)
    vol_bedsusT=sum(vol_bedsus)
    vol_zirT=sum(vol_zir)
    
    #                           Check for mass balance                            #
    
#    check=len(x)-sum(npixels)
#    if abs(check)>0:
#        print('Calculation problem, ' +str(check)+ ' pixels lost')
#    check=vol_bedsusT-len(x)*eros*pixel_size*pixel_size
#    if abs(check)>0:
#        print('Mass balance problem, ' +str(check)+ ' m3 lost')
#    
    #                            Calculate bedload ratio                          #
    
    bedloadratio = vol_bedT/vol_bedsusT*100
#    print('At the outlet considered, bedload represents ' \
#    + str("%.1f" % bedloadratio) + ' % of the sediment.')
#    
    #              PERCENT arrays that contain data for each litho                #
    
    PC_pixels=[]; PC_bed=[]; PC_sus=[]; PC_bedsus=[]; PC_zir=[]; 
    for i in range(len(npixels)):
        PC_pixels.append(100*float(npixels[i])/float(len(x)))
        PC_bed.append(100*vol_bed[i]/vol_bedT)
        PC_sus.append(100*vol_sus[i]/vol_susT)
        PC_bedsus.append(100*vol_bedsus[i]/vol_bedsusT)
        PC_zir.append(100*vol_zir[i]/vol_zirT)    
            
    #                           Displaying the results                            #
            
#    for i in range(len(npixels)):
#        print('LITHOLOGY '+str(i+1)+':')
#        print('- makes up '+str("%.1f" % PC_pixels[i])+'% of catchment.')
#        print('- makes up '+str("%.1f" % PC_bed[i])+'% of bedload, '\
#        +str("%.1f" % PC_sus[i])+'% of suspended load and '\
#        +str("%.1f" % PC_bedsus[i])+'% of total load.')
#        print( '- makes up '+str("%.1f" % PC_zir[i])+'% of zircons found in sands.')
    
    phig = np.array([15.5, 10, 0, 13.5, 61])
    #np.array([32, 12, 0, 3, 53])
    print(erod)
    cost = 100*np.sum( (phig - np.array(PC_zir))**2 )/(5*np.sum(phig**2))
    print(PC_zir)
    print(cost)
    print()
    return cost
    #return PC_zir[:]


#PDF['Art31'] = 0.464*PDF['A'] + 0.128*PDF['C'] + 0.173*PDF['F'] + 0.235*PDF['H'] # Max grain-size (0.35 for TTS and 0.15 for the rest)
#phi_ami2 = [0.654, 0.083, 0.111, 0.152] # Max erosion (5.1 for TTS and 1 for the rest)
#PDF['Art31'] = 0.565*PDF['A'] + 0.104*PDF['C'] + 0.14*PDF['F'] + 0.191*PDF['H'] # Max grain-size (0.5 for TTS and 0.15 for the rest)
#phi_abr = [0.525, 0.129, 0.158, 0.188] # Max abrasion (31 for TTS and 0.15 for the rest)
#PDF['Art'] = 0.492*PDF['A'] + 0.122*PDF['C'] + 0.164*PDF['F'] + 0.222*PDF['H'] # Max fertility (8.1 for TTS and 3.1 for rest)
#PDF['Art'] = 0.155*PDF['A'] + 0.10*PDF['C'] + 0.135*PDF['F'] + 0.61*PDF['H'] # Max grain-size (0.5 for LH and 0.15 for the rest)

from scipy.optimize import minimize
initial_phi = np.ones(5)*2#[0.1, 0.1, 0.1, 0.1, 10]
bounds = tuple([(0.15, 31)]*len(initial_phi))
#constraints = ({'type':'eq', 'fun': lambda x: sum(x)-1},)
res = minimize(cost_function_mistmatch, initial_phi, method='SLSQP',#constraints=constraints,
               bounds=bounds, tol=1e-3, options={'disp':True})
lowest_cost = res.fun
best_phi = res.x

    
from scipy.optimize import minimize
initial_erod = np.ones(5)*5 #[0.1, 0.1, 0.1, 0.1, 10]

bounds = tuple([(0.15,31)]*len(initial_erod))
lowest_cost = 100
for i in range(10):
    initial_erod = np.random.uniform(0.15,31, size=(5,1))
    #constraints = ({'type':'eq', 'fun': lambda x: sum(x)-1},)
    res = minimize(cost_function_mistmatch, initial_erod, method='SLSQP',#constraints=constraints,
               bounds=bounds, tol=1e-3, options={'disp':True})#, 'maxiter':20
    if res.fun < lowest_cost:
        lowest_cost = res.fun
        best_phi = res.x
        print('Found new min at %s with %.2f%% fit' % (best_phi, lowest_cost))
with open('abrasion_test_vs_abrasion.txt', 'w') as f:
    f.write("Abrasion reproducing abrasion: %s\n%.2f\n%s" % (best_phi, lowest_cost, PC_zir))

initial_phi = np.ones(5)*2#[0.1, 0.1, 0.1, 0.1, 10]
bounds = tuple([(0.1,31)]*len(initial_phi))
lowest_cost = 100
for i in range(10):
    initial_phi = np.random.uniform(0.15,31, size=(5,1))
    res = minimize(cost_function_mistmatch, initial_phi, method='SLSQP',#constraints=constraints,
                   bounds=bounds, tol=1e-3, options={'disp':True})
    if res.fun < lowest_cost:
        lowest_cost = res.fun
        best_phi = res.x
        print('Found new min at %s with %.2f%% fit' % (best_phi, lowest_cost))