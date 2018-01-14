# -*- coding: utf-8 -*-

"""
 
Pebble abrasion model (pABRASIONmodel.py) calculates mixing proportions of
source material in modern river sediments.

 Description

   pABRASIONmodel.py is a Python 3 module contaning functions ('analysis', 
   'cost_function_LSQ' and 'minimisation') that calculate area, bedload, 
   suspended load, and zircon mixing proportion at the river outlet. 
   It also returns how well a factor being tested is able to isolatedly 
   reproduce the distortion caused by another factor (target) on the zircon 
   mixing proportions in modern river sediments.

   For a detailed explanation of the equations, please refer to the work of
   Lavarini et al. [2018a].

 Reference:
     
     Lavarini et al. [2018]. Does pebble abrasion influence detrital age 
     population statistics? A numerical investigation of natural data sets. 
     Journal of Geophysical Research - Earth Surface.

 Authors: C. Lavarini (lavarini.c@gmail.com), C. A. da Costa Filho & M. Attal
 The University of Edinburgh, UK
 December of 2017
      
 """
## Example of inputs for the function 'analysis'. 
# Working directory:
DataDirectory = r'D:/Chrystiann Lavarini/Universidade/Doutorado/PhD Project/MODEL/' 

# Watershed files extracted through GIS:
FileName = 'MARS_90_DEM_PT.dbf' # Elevation.
FileName1 = 'F_ACC_MARS_90_PT.dbf' # Flow accumulation.
FileName2 = 'EXT_GEO_MARS_K_BASIN.dbf'  # Lithological unity.
FileName3 = 'F_LEN_MARS_90_PT.dbf' # Flow length.

pixel_size=90 # Pixel size

eros = [1, 1, 1, 1, 1]   # Erosion rates (mm/year)
erod = [0, 0, 0, 0, 0]   # Abrasion rates (mass loss/km)
fert = [1, 1, 1, 1, 1]   # Fertility (relative concentration)
source_1mm = [0.25, 0.25, 0.25, 0.25, 0.25]  # Gravel supply (%)

# Source units:
litho1 = 'Tethyan series'
litho2 = 'Formation II/III'
litho3 = 'Manaslu granite'
litho4 = 'Formation I'
litho5 = 'Lesser Himalaya'

def analysis(DataDirectory, FileName, FileName1, FileName2, FileName3, \
                pixel_size, eros, erod, fert, source_1mm, litho1, litho2,  \
                litho3, litho4, litho5):
    """
    
    This function executes a code that calculates area, 
    bedload, suspended load and zircon mixing proportions at the river's 
    outlet. 
    
    Exposure area, lithological units, elevation, flow length, flow 
    accumulation, abrasion rate, gravel supply, mineral concentration, and 
    erosion rate are user-dependent arguments. It returns area, bedload, 
    suspended load and percentage of zircon from the source units. 

 Input arguments

   DataDirectory = Working directory. Type: String
   FileName = DEM (contaning x, y and z). Type: *.dbf
   FileName2 = Flow accumulation (contaning x, y). Type: *.dbf
   FileName3 = Lithological units (contaning x, y). Type: *.dbf
   FileName4 = Flow length (contaning x, y). Type: *.dbf
   pixel_size = DEM pixel size (in square meters). Type: Float
   eros = Erosion rates (in mm/year). Type: List
   erod = Abrasion rates (in mass loss/km). Type: List
   fert = Fertility (in relative concentration). Type: List
   source_1mm = Sand supply (in % 10^(-2)). Type: List
   
   For a detailed verification of input data, please, see the input arguments
   available in the 'Examples' folder at GitHub.

 Output arguments 

   'pABRASIONmodel_outputs.txt' and 'pABRASION_output.xlsx' contain area, 
   percentage of bedload, percentage of suspended load and percentage of 
   zircon from each lithological unit in the mixed river sediment.

        
 Example


    FileName = 'DEM.dbf' 
    FileName1 = 'F_ACC.dbf' 
    FileName2 = 'GEOLOGY.dbf' 
    FileName3 = 'F_LEN.dbf'                                                                  
    pixel_size = 90
    eros = [1, 2, 5.1, 0.1, 3] 
    erod = [10, 2, 15, 0.1, 30] 
    fert = [1, 2, 1, 0.1, 3]
    source_1mm = [0.25, 0.25, 0.25, 0.25, 0.25]
    
    The lists 'eros', erod, fert and source_1mm are ordered. The first element
    refers to lithology 1, the second to lithology 2, etc.
    
   For a detailed explanation of the equations, please refer to the work of
   Lavarini et al. [2018a].
    
    """
    # Import modules required to perform calculation:
    import numpy as np
    from dbfread import DBF
    import pandas as pd

    # Creates series of arrays to perform calculation:                                                                                                   
    PC_zir_ALL = [] 
    bedloadratio_ALL = [] 
    PC_sus_ALL = [] 
    PC_bed_ALL = [] 
    PC_pixels_ALL = [] 
       
    # Coordinates (x, y) and elevation (z) from FileName:    
    dbffile = DBF(DataDirectory + FileName)
    global x

    x = []
    y = []
    z = []
    
    for record in dbffile:
        x.append(record['x'])
        y.append(record['y'])
        z.append(record['grid_code'])
    
    # Area of every lithological unit from FileName1:  
    dbffile = DBF(DataDirectory + FileName1)
    
    global area
    area = []
    
    for record in dbffile:
        area.append(record['grid_code'])
    
    # The code of every lithological unit from FileName:              
    dbffile = DBF(DataDirectory + FileName2)
    global litho
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
          
    # Flow length is collected (from FileName3):                        
    dbffile = DBF(DataDirectory + FileName3)
    global length
    length = []
    
    for record in dbffile:
        length.append(record['grid_code'])    
       
    #Set parameters for model: 
    global n_lithos
    n_lithos = len(set(LITHO))
    
    litho_names = []
    
    for i in LITHO:
        if LITHO[i] ==1:
            litho_names.append(str(litho1)) # Give the name of this formation.
        elif LITHO[i] ==2:
            litho_names.append(str(litho2)) # Give the name of this formation.
        elif LITHO[i] ==3:
            litho_names.append(str(litho3)) # Give the name of this formation.
        elif LITHO[i] ==4:
            litho_names.append(str(litho4)) # Give the name of this formation.
        else:
            litho_names.append(str(litho5)) # Give the name of this formation.        
    
    max_area=0.0; outletID=0;
    n_lines = len(x)   # get the number of pixels (=number of data)             
    for i in range(0, n_lines):
        if area[i] > max_area:  #  finding the outlet node and basin total area 
            max_area=area[i]; outletID=i 
            
    # Printing and iterating:                           
    print() 
    print( "River outlet")  
    print()            
    print( "Outlet node ID is "+str(outletID)+" at x="+str(x[outletID])+" m and \
    y="+str(y[outletID])+" m; drainage area is "+str(area[outletID])+" pixels.")
    
    # Setting the length of all nodes to outlet node:    
    length2outlet=[]
    for i in range(len(length)):
        length2outlet.append(length[i]-length[outletID])
        # Checking for errors:
        if length[i]-length[outletID]<0.0:
          print( "Error: all nodes should be upstream of outlet") 
    
    # Lists to store data                                                                      
    npixels=[]; vol_bed=[]; vol_sus=[]; vol_bedsus=[]; vol_zir=[]; 
    vol_bedT=0; vol_susT=0; vol_bedsusT=0; vol_zirT=0; 
    
    for i in range(1, n_lithos+1):    # proceed for each litho
        count=0; bedload=0; suspended=0; 
        print( 'Proceeding with litho number '+str(i))
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
           
        npixels.append(count)   # arrays contains number of pixels for each litho 
        vol_bed.append(bedload) # volume of bedload for each litho (m3/yr)        #
        vol_sus.append(suspended) # volume of suspended load for litho (m3/yr)    #
        total_volumeforlitho=bedload+suspended # contains total volume of sediment#   
        vol_bedsus.append(total_volumeforlitho)  # for each litho (m3/yr)         #
        total_zirconforlitho=suspended*fert[i-1] # contains volume of zircon in   #
        vol_zir.append(total_zirconforlitho)     # sand for each litho (m3/yr)    # 
       
    #                           Check for mass balance                            #
        check=suspended+bedload-(count*eros[i-1]*pixel_size*pixel_size)
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
    check=vol_bedsusT-len(x)*eros[i-1]*pixel_size*pixel_size
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
    
    ## Saving a text file with the results:
    
    with open('pABRASIONmodel_outputs', 'w') as f:
        for i in range(len(npixels)):
            f.write('LITHOLOGY '+str(i+1)+':\n')
            f.write('- makes up '+str("%.2f" % PC_pixels[i])+'% of catchment.\n')
            f.write('- makes up '+str("%.2f" % PC_bed[i])+'% of bedload, '\
                    +str("%.2f" % PC_sus[i])+'% of suspended load and '\
                    +str("%.2f" % PC_bedsus[i])+'% of total load.\n')
            f.write( '- makes up '+str("%.2f" % PC_zir[i])+'% of zircons found in sands.\n')
            f.write( '-------------------------------------------------------------------\n')
            f.write( ' The input arguments used in this experiment for lithology '+str(i+1)+ ' were: \n')
            f.write( 'Abrasion rates: ' + str("%.2f" % erod[i]) + '\n')       
            f.write( 'Erosion rates: ' + str("%.2f" % eros[i]) + '\n')
            f.write( 'Specific mineral fertility: ' + str("%.2f" % fert[i]) + '\n')
            f.write( 'Gravel supply: ' + str("%.2f" % source_1mm[i]) + '\n')
            
    ## Saving a table file with the results:
    
    # Create columns for the Excel:
    l1 = ['Area [% of catchment]', 'Bedload [%]', 'Suspended load [%]', \
          'Total load', 'Zircon [%]', 'Abrasion rates [% mass loss/km]', \
          'Erosion rates [mm/year]', 'Zircon fertility [grain/g]', 'Gravel supply [%]']
    l2 = []
    l3 = []
    l4 = []
    l5 = []
    l6 = []
    
    # Iterate and store the inputs and outputs:
    for i in range(len(npixels)):
        if i == 0:        
            l2.append(PC_pixels[i]), l2.append(PC_bed[i]), l2.append(PC_sus[i]), \
            l2.append(PC_bedsus[i]), l2.append(PC_zir[i]), l2.append(erod[i]), \
            l2.append(eros[i]), l2.append(fert[i]), l2.append(source_1mm[i])
    
        elif i == 1:        
            l3.append(PC_pixels[i]), l3.append(PC_bed[i]), l3.append(PC_sus[i]), \
            l3.append(PC_bedsus[i]), l3.append(PC_zir[i]), l3.append(erod[i]), \
            l3.append(eros[i]), l3.append(fert[i]), l3.append(source_1mm[i])
    
        elif i == 2:        
            l4.append(PC_pixels[i]), l4.append(PC_bed[i]), l4.append(PC_sus[i]), \
            l4.append(PC_bedsus[i]), l4.append(PC_zir[i]), l4.append(erod[i]), \
            l4.append(eros[i]), l4.append(fert[i]), l4.append(source_1mm[i])
    
        elif i == 3:        
            l5.append(PC_pixels[i]), l5.append(PC_bed[i]), l5.append(PC_sus[i]), \
            l5.append(PC_bedsus[i]), l5.append(PC_zir[i]), l5.append(erod[i]), \
            l5.append(eros[i]), l5.append(fert[i]), l5.append(source_1mm[i])
    
        elif i == 4:        
            l6.append(PC_pixels[i]), l6.append(PC_bed[i]), l6.append(PC_sus[i]), \
            l6.append(PC_bedsus[i]), l6.append(PC_zir[i]), l6.append(erod[i]), \
            l6.append(eros[i]), l6.append(fert[i]), l6.append(source_1mm[i])        
    
    # Create a data frame to store the output:
    df = pd.DataFrame({'Variable': l1, litho1: l2, litho2: l3, litho3: l4, \
                    litho4: l5, litho5: l6})
    # Create a Pandas Excel writer using XlsxWriter as the engine:
    writer = pd.ExcelWriter('pABRASION_output00.xlsx', engine='xlsxwriter')
    
    # Convert the dataframe to an XlsxWriter Excel object:
    df.to_excel(writer, sheet_name='Sheet1', index=False)
    
    # Close the Pandas Excel writer and output the Excel file:
    writer.save()

    return print('Your experiment has finished. Please, check your output ' \
                 + 'files in the working directory.')

## Cost function for minimisation:

def cost_function_LSQ(erod):
    
    """
 
    The cost function of least-squares ('cost_function_LSQ') estimates how 
    well a factor being tested is able to isolatedly reproduce the distortion 
    caused by another factor (target) on the zircon mixing proportions.

     Description
    
       'cost_function_LSQ'is based on the characteristics of sediment source 
       units. Abrasion rate, gravel supply, mineral concentration, erosion rate 
       and target zircon mixing proportion are user-dependent arguments. It 
       returns the zircon mixing proportions made by the tested factors and 
       how well it replicates the target factors (i.e., a cost function).
       
       For a detailed explanation of the equations, please refer to the work of
       Lavarini et al. [2018a].
    
     Input arguments
    
       eros = Erosion rates (in mm/year). Type: List
       erod = Abrasion rates (in mass loss/km). Type: List
       fert = Fertility (in relative concentration). Type: List
       source_1mm = Sand supply (in % 10^(-2)). Type: List
       
       Caution note: This function needs to be altered when the user changes
       the testing factors. Moreover, before a factor being tested, it must 
       be removed from the code, otherwise it overwrites and the minimisation
       fails (Please, see the topic 'Example' below). 
           
     Output arguments 
    
       It returns the cost function that can be further minimised.
    
            
     Example
    
        Target factor: erod (i.e., abrasion rate)
        
        # Because the target factor is abrasion (erod), the list of abrasion 
        rates in the function (erod) needs to be removed (or commented). 
        Otherwise, it generates error in the minimisation.       
        
        # A comment is made for erod:
        #erod = [2, 0, 0, 0, 0]   # Abrasion rates (mass loss/km)
        eros = [1, 1, 1, 1, 1]   # Erosion rates (mm/year)
        fert = [1, 1, 0, 1, 1]   # Fertility (relative concentration)
        source_1mm = [0.25, 0.25, 0.25, 0.25, 0.25]  # Gravel supply (%)

        # The values in the target zircon mixing proportions (phig) are assigned:        
        phig = np.array([52.5, 12.9, 0, 15.8, 18.8]) 
    
        # The testing variable needs to be assigned here too for displaying:
        print(erod) 
    
    
     Reference:
         
         Lavarini et al. [2018]. Does pebble abrasion influence detrital age 
         population statistics? A numerical investigation of natural data sets. 
         Journal of Geophysical Research - Earth Surface.
    
     Authors: C. Lavarini (lavarini.c@gmail.com), C. A. da Costa Filho & M. Attal
     The University of Edinburgh, UK
     December of 2017
      
 """
    # A "gentle" remind: 
    print('Be aware of the caution note in the docstring of this function!')
    import numpy as np

    ## The inputs below must be changed/removed accordingly when needed.    
    eros = [1, 1, 1, 1, 1]   # Erosion rates (mm/year)
    #erod = [2, 0, 0, 0, 0]   # Abrasion rates (mass loss/km)
    fert = [1, 1, 0, 1, 1]   # Fertility (relative concentration)
    source_1mm = [0.25, 0.25, 0.25, 0.25, 0.25]  # Gravel supply (%)

    ## The code between the lines 428 and 470 are the same as the function 'analysis'.    
    max_area=0.0; outletID=0;
    n_lines = len(x)   
    for i in range(0, n_lines):
        if area[i] > max_area: 
            max_area=area[i]; outletID=i 
    
    length2outlet=[]
    for i in range(len(length)):
        length2outlet.append(length[i]-length[outletID])
        if length[i]-length[outletID]<0.0:
          print( "Error: all nodes should be upstream of outlet") 
                                                                          
    npixels=[]; vol_bed=[]; vol_sus=[]; vol_bedsus=[]; vol_zir=[]; 
    vol_bedT=0; vol_susT=0; vol_bedsusT=0; vol_zirT=0; 
    
    for i in range(1, n_lithos+1): 
        count=0; bedload=0; suspended=0; 

        for j in range(len(x)):     
            if litho[j]==i:       
                bedload=bedload + (1-source_1mm[i-1])*eros[i-1]*pixel_size*pixel_size*\
                    np.exp(-length2outlet[j]/1000*erod[i-1]/100)
                suspended=suspended + source_1mm[i-1]*eros[i-1]*pixel_size*pixel_size +\
                    (1-source_1mm[i-1])*eros[i-1]*pixel_size*pixel_size*\
                    (1-np.exp(-length2outlet[j]/1000*erod[i-1]/100))
                count = count + 1
            
        npixels.append(count)   
        vol_bed.append(bedload) 
        vol_sus.append(suspended) 
        total_volumeforlitho=bedload+suspended 
        vol_bedsus.append(total_volumeforlitho)  
        total_zirconforlitho=suspended*fert[i-1] 
        vol_zir.append(total_zirconforlitho)    
       
        check=suspended+bedload-(count*eros[i-1]*pixel_size*pixel_size)
    
    vol_bedT=sum(vol_bed)
    vol_susT=sum(vol_sus)
    vol_bedsusT=sum(vol_bedsus)
    vol_zirT=sum(vol_zir)
        
    bedloadratio = vol_bedT/vol_bedsusT*100
        
    # The zircon mixing proportion (PC_zir) is assigned as global
    global PC_zir
    
    ##  PERCENT arrays that contain data for each lithology   
    PC_pixels=[]; PC_bed=[]; PC_sus=[]; PC_bedsus=[]; PC_zir=[]; 
    for i in range(len(npixels)):
        PC_pixels.append(100*float(npixels[i])/float(len(x)))
        PC_bed.append(100*vol_bed[i]/vol_bedT)
        PC_sus.append(100*vol_sus[i]/vol_susT)
        PC_bedsus.append(100*vol_bedsus[i]/vol_bedsusT)
        PC_zir.append(100*vol_zir[i]/vol_zirT)    
    
    ### The Target mixing proportion and testing variable must be changed if needed.  
    ## Target mixing proportion
    global phig      
    phig = np.array([52.5, 12.9, 0, 15.8, 18.8]) 
    
    # Display the testing variable
    print(erod) 
    
    ## Cost function    
    cost = 100*np.sum( (phig - np.array(PC_zir))**2 )/(5*np.sum(phig**2))

    ## Displaying results        
    print(PC_zir)
    print(cost)
    print()
    return cost


# Import minimisation function
def minimisation(cost_function_LSQ):

    """
 
    The minimisation function ('minimisation') estimates how well a 
    factor being tested is able to isolatedly reproduce the distortion caused 
    by another factor (target) on the zircon mixing proportions.

     Description
    
       The minimisation function takes the 'cost_function_LSQ' as an argument.
       It returns the zircon mixing proportions made by the tested factors and 
       how well it replicates the target factors (i.e., a cost function).
       
       For a detailed explanation of the equations, please refer to the work of
       Lavarini et al. [2018a].
    
     Input arguments
         
        'cost_function_LSQ'= Cost function of least-squares. Type: Function
                  
     Output arguments 
    
        It returns a file 'zMIXTUREminimisation.xlsx'and with the zircon mixing 
        proportions made by the tested factors from the function 
        'cost_function_LSQ' and how well it replicates the target factors 
        (i.e., a cost function itself).
    
     Reference:
         
         Lavarini et al. [2018]. Does pebble abrasion influence detrital age 
         population statistics? A numerical investigation of natural data sets. 
         Journal of Geophysical Research - Earth Surface.
    
     Authors: C. Lavarini (lavarini.c@gmail.com), C. A. da Costa Filho & M. Attal
     The University of Edinburgh, UK
     December of 2017
      
 """
    # Import the minimize function from SciPy
    from scipy.optimize import minimize
    import numpy as np
    import pandas as pd
    
    # Initial values to be used in the minimisation
    initial_phi = np.ones(5)*2
    
    # Boundary conditions 
    bounds = tuple([(0.15, 31)]*len(initial_phi))
    
    # Minimisation function
    res = minimize(cost_function_LSQ, initial_phi, method='SLSQP',
                   bounds=bounds, tol=1e-3, options={'disp':True})
    lowest_cost = res.fun
    best_phi = res.x
    
    # Export the results
    with open('zMIXTUREminimisation', 'w') as f:
        f.write("Minimisation result: %s\n%.2f\n%s" % (best_phi, lowest_cost, PC_zir))
        
        ## Saving a table file with the results:
    
    # Create columns for the Excel:
    l1 = ['Testing variable [best-fit]', 'Zircon proportion [%, made by the testing variable]', \
          'Target zircon proportion [%]', 'Cost function [least-square]']
    l2 = []
    l3 = []
    l4 = []
    l5 = []
    l6 = []
    
    # Iterate and store the inputs and outputs:
    for i in range(len(erod)):
        if i == 0:        
            l2.append(best_phi[i]), l2.append(PC_zir[i]), l2.append(phig[i]), \
            l2.append(lowest_cost)
        elif i == 1:        
            l3.append(best_phi[i]), l3.append(PC_zir[i]), l3.append(phig[i]), \
            l3.append(lowest_cost)
    
        elif i == 2:        
            l4.append(best_phi[i]), l4.append(PC_zir[i]), l4.append(phig[i]), \
            l4.append(lowest_cost)
    
        elif i == 3:        
            l5.append(best_phi[i]), l5.append(PC_zir[i]), l5.append(phig[i]), \
            l5.append(lowest_cost)
    
        elif i == 4:        
            l6.append(best_phi[i]), l6.append(PC_zir[i]), l6.append(phig[i]), \
            l6.append(lowest_cost)    
    
    # Create a data frame to store the output:
    df = pd.DataFrame({'Variable': l1, litho1: l2, litho2: l3, litho3: l4, \
                    litho4: l5, litho5: l6})
    # Create a Pandas Excel writer using XlsxWriter as the engine:
    writer = pd.ExcelWriter('zMIXTUREminimisation.xlsx', engine='xlsxwriter')
    
    # Convert the dataframe to an XlsxWriter Excel object:
    df.to_excel(writer, sheet_name='Sheet1', index=False)
    
    # Close the Pandas Excel writer and output the Excel file:
    writer.save()

    return print('Your experiment has finished. Please, check your output ' \
                 + 'files in the working directory.')

