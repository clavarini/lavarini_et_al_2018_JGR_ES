# -*- coding: utf-8 -*-

"""
 
Zircon ages module (zAGES.py) generates probability density functions (PDFs) 
and provides tools for analysing them.

 Description

   zAGES.py is a Python 3 module contaning functions ('generate_pdf', 
   'pdfMISMATCH',  'zMIX', 'zSTATS', 'qqPLOT', and 'pdfCROSS') that generates
   probability density functions (PDFs) from U-Pb ages. It also has tools for 
   discovering the best-fit of zircon mixing based on the U-Pb ages as well as
   comparing the age samples (area mismatch, similarity, likeliness, 
   Kolmogorov-Smirnov, Q-Q plot and PDF crossplots). 

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
 
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import numpy as np
import os
from numpy import pi, exp, sqrt
from scipy import stats
from scipy.optimize import minimize
from scipy.integrate import cumtrapz
import statsmodels.api as sm
blu, grn, red, prp, yel, cyn = sns.color_palette()

sns.set(style='ticks', font_scale=1.5, 
        rc={'lines.linewidth': 2.5,
            'figure.figsize': (10, 8),
            'text.usetex': False,
            # 'font.family': 'sans-serif',
            # 'font.sans-serif': 'Optima LT Std',
        })

figpath = os.path.join('K_AGES/', 'marsyandi_')

plt.ion()
sheets = ['A', 'C', 'F', 'H', 'K']
filename = 'MARSYANDI_SYNTHETIC_1.xlsx'
objects = ('A', 'C', 'K', 'H', 'F')
objects2 = ('A', 'C', 'H', 'F')

PDF = {}
SIZES = {}

def generate_pdf(filename, sheet, smooth=80, plot=False):
    
    """
    Generate PDFs  ('generate_pdf') is a function that generates probability 
    density plots (PDPs) and bar plots of U-Pb ages.
    
     Description
    
       'generate_pdf' is a function that reads Excel spreadsheets with U-Pb
       ages and then creates a histogram and a probability density plot (PDP)
       for each U-Pb age list.
       
       For a detailed explanation of the equations, please refer to the work of
       Lavarini et al. [2018a].

     Input arguments
    
       figpath = The path to save the plots. Type: String
       filename =. The name of the Excel spreadsheet. Type: String
       sheet = The name of the sheet under analysis. Type: List
       smooth = Smoothing window interval. Type: Integer
       plot = True or False. Type: Boolean
                  
     Output arguments 
    
       It returns the histogram and PDP of U-Pb ages.
             
     Example
    
        figpath = os.path.join('K_AGES/', 'marsyandi_')
        
        plt.ion()
        sheets = ['A', 'C', 'F', 'H', 'K']
        filename = 'MARSYANDI_SYNTHETIC_1.xlsx'
        objects = ('A', 'C', 'K', 'H', 'F')
        objects2 = ('A', 'C', 'H', 'F')
        
        PDF = {}
        SIZES = {}
        
        for sheet in sheets:
            PDF[sheet], SIZES[sheet] = generate_pdf(filename, \
            sheet, smooth=80, plot=True)

     Reference:
         
         Lavarini et al. [2018]. Does pebble abrasion influence detrital age 
         population statistics? A numerical investigation of natural data sets. 
         Journal of Geophysical Research - Earth Surface.
    
     Authors: C. Lavarini (lavarini.c@gmail.com), C. A. da Costa Filho & M. Attal
     The University of Edinburgh, UK
     December of 2017
    
    """
    df = pd.read_excel(filename, sheetname = sheet, skiprows=(0,))

    # Create single-column with the ages
    K = df[['Age (Ma)', 'Error (Ma)']]
    K.columns = ['m', 's']
    K.name = 'Sample '+sheet

    # Create a column for Age (Ma) and one for Error (Ma)
    K['m'] = df['Age (Ma)'] 
    K['s']  = df['Error (Ma)'] 

    # Plot U-Pb histogram
    if plot:
        f, ax = plt.subplots()
        plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        K['m'].hist(bins=50, ax=ax, color='0.6')
        ax.set_xlim([0, 4000])
        ax.set_xlabel('Synthetic U-Pb age [Ma]')
        ax.set_ylabel('Frequency')
        
        # Save plots in the chosen path.
        f.savefig(figpath+'histo_'+sheet+'.pdf', bbox_inches='tight')
        f.savefig(figpath+'histo_'+sheet+'.png', bbox_inches='tight')
    sample_size = len(K['m'])
    
    # Calculate probability density function (PDF)
    x = np.linspace(0, 4000, 4001) # Ma
    pdf = np.zeros(len(x))
    for i, row in K.iterrows():
        # Create PDF for each sample
        mu_i, sigma_i = row['m'], row['s']
        pdf += exp( -(x - mu_i)**2/(2*sigma_i**2) )/(2*sigma_i*sqrt(2*pi))
    
    # Plot PDF for each sample
    if plot:
        f, ax = plt.subplots()
        plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        ax.fill_between(x, 0, pdf, color='k')
        ax.set_xlabel('Synthetic U-Pb age [Ma]')
        ax.set_ylabel('Relative probability')
        
        # Save plots in the chosen path.
        f.savefig(figpath+'pdf_spiky_'+sheet+'.pdf', bbox_inches='tight')
        f.savefig(figpath+'pdf_spiky_'+sheet+'.png', bbox_inches='tight')

    # Create DataFrame for PDF to aid smoothing
    pdf = pd.DataFrame(pdf, index=x)
    pdf = pdf.rolling(window=smooth, center=True, win_type='boxcar').mean().fillna(0)
    pdf /= pdf.sum()
    if plot:
        f, ax = plt.subplots()
        plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        ax.fill_between(x, 0, pdf.as_matrix().flatten(), color='k')
        ax.set_xlabel('Synthetic U-Pb age [Ma]')
        ax.set_ylabel('Relative probability')
        textstr = sheet
        props = dict(boxstyle='round', facecolor='white', alpha=1.0)
        ax.text(0.6, 0.5, textstr, transform=ax.transAxes,
        verticalalignment='center', bbox=props)
        
        # Save plots in the chosen path.
        f.savefig(figpath+'pdf_smooth_'+sheet+'.pdf', bbox_inches='tight')
        f.savefig(figpath+'pdf_smooth_'+sheet+'.png', bbox_inches='tight')
    return pdf, sample_size

for sheet in sheets:
    PDF[sheet], SIZES[sheet] = generate_pdf(filename, sheet, smooth=80, plot=False)

def pdfMISMATCH(phi):
    
    """
 
    The area mismatch ('pdfMISMATCH') function estimates the percentage of areal
    mismatch between a PDF made of mixture of zircons of single sediment sources
    and a downstream mixed sample.

     Description
    
       'pdfMISMATCH'is a cost function for finding the best mixing proportion of 
       zircons between single sediment sources and a downstream mixed sample.
       It is based on the area mismatch between an U-Pb artificial PDF and a 
       U-Pb real mixed downstream PDF sample, as suggested by Amidon et al. 
       [2005a].

       For a detailed explanation of the equations, please refer to the work of
       Lavarini et al. [2018a] or Amidon et al. [2005a].
    
     Input arguments
    
       phi = Initial values to be used as mixing proportion. Type: List
           
     Output arguments 
    
       It returns the area mismatch between the artificial and real PDFs
       that can be further minimised.
       
       There is no file exported.
    
     Reference:
         
         Lavarini et al. [2018]. Does pebble abrasion influence detrital age 
         population statistics? A numerical investigation of natural data sets. 
         Journal of Geophysical Research - Earth Surface.
         
         Amidon et al. [2005a]. Construction of detrital mineral populations: 
             Insights from mixing of U-Pb zircon ages in Himalayan rivers.
             Basin Research.
    
     Authors: C. Lavarini (lavarini.c@gmail.com), C. A. da Costa Filho & M. Attal
     The University of Edinburgh, UK
     December of 2017
      
     """
    f = phi[0]*PDF['A'] + phi[1]*PDF['C'] + \
    phi[2]*PDF['F'] + phi[3]*PDF['H']
    x = f.index
    cost = np.trapz(abs(f-PDF['K']).as_matrix().flatten(),
              dx=x[2]-x[1])/2
#    print('The area mismatch [%] is: ', cost)
    return cost*100

def zMIX(pdfMISMATCH):

    """
 
    The zircon mixing proporiton ('zMIX') function estimates the best-fit 
    mixture of zircons between single sediment sources and a downstream mixed 
    sample.

     Description
    
       'zMIX' uses a cost function (area mismatch) for finding the best 
       mixing proportion of zircons between single sediment sources and a 
       downstream mixed sample.
       It is based on the area mismatch between an U-Pb artificial PDF and a 
       U-Pb real mixed downstream PDF sample, as suggested by Amidon et al. 
       [2005a].

       For a detailed explanation of the equations, please refer to the work of
       Lavarini et al. [2018a] or Amidon et al. [2005a].
    
     Input arguments
    
       pdfMISMATCH = It is the area mismatch function. Type: Function
           
     Output arguments 
    
       It returns an Excel file with the best-fit zircon mixing proportion and
       the area mismatch.
    
     Reference:
         
         Lavarini et al. [2018]. Does pebble abrasion influence detrital age 
         population statistics? A numerical investigation of natural data sets. 
         Journal of Geophysical Research - Earth Surface.
         
         Amidon et al. [2005a]. Construction of detrital mineral populations: 
             Insights from mixing of U-Pb zircon ages in Himalayan rivers.
             Basin Research.
    
     Authors: C. Lavarini (lavarini.c@gmail.com), C. A. da Costa Filho & M. Attal
     The University of Edinburgh, UK
     December of 2017
      
     """
     
    initial_phi = [0, 1, 0, 0]
    bounds = tuple([(0,1)]*len(initial_phi))
    constraints = ({'type':'eq', 'fun': lambda x: sum(x)-1},)
    res = minimize(pdfMISMATCH, initial_phi, method='SLSQP', bounds=bounds, constraints=constraints, tol=1e-10)
    lowest_cost = res.fun
    global best_phi
    best_phi = res.x
    for i in range(30):
        initial_phi = np.random.uniform(0, 1, size=(len(initial_phi),))
        initial_phi /= sum(initial_phi)
        bounds = tuple([(0,1)]*len(initial_phi))
        res = minimize(pdfMISMATCH, initial_phi, method='SLSQP', bounds=bounds, constraints=constraints, tol=1e-10)
        if res.fun < lowest_cost:
            lowest_cost = res.fun
            best_phi = res.x
            print()
            print('Area mismatch [%]: ', lowest_cost)
            print('Zircon mixing proportion [% $10^2$]:  ', best_phi)
    
            ## Saving a table file with the results:
    
    # Create columns for the Excel:
    l1 = ['Zircon proportion [% 10^2]', \
           'Area mismatch [%]']
    l2 = []
    l3 = []
    l4 = []
    l5 = []
    
    # Iterate and store the inputs and outputs:
    for i in range(len(initial_phi)):
        if i == 0:        
            l2.append(best_phi[i]),  \
            l2.append(lowest_cost)
        elif i == 1:        
            l3.append(best_phi[i]), \
            l3.append(lowest_cost)
    
        elif i == 2:        
            l4.append(best_phi[i]),  \
            l4.append(lowest_cost)
    
        elif i == 3:        
            l5.append(best_phi[i]),  \
            l5.append(lowest_cost)
    
    
    # Create a data frame to store the output:
    df = pd.DataFrame({'Variable': l1, 'litho1': l2, 'litho2': l3, 'litho3': l4, \
                    'litho4': l5})
    # Create a Pandas Excel writer using XlsxWriter as the engine:
    writer = pd.ExcelWriter('zMIXTUREages.xlsx', engine='xlsxwriter')
    
    # Convert the dataframe to an XlsxWriter Excel object:
    df.to_excel(writer, sheet_name='Sheet1', index=False)
    
    # Close the Pandas Excel writer and output the Excel file:
    writer.save()

    return print('Your experiment has finished. Please, check your output ' \
                 + 'files in the working directory.')
    
Ami = best_phi[0]*PDF['A'] + best_phi[1]*PDF['C'] + \
         best_phi[2]*PDF['F'] + best_phi[3]*PDF['H']
Lav = PDF['K']

def zSTATS(Ami, Lav):
    
    """
    zSTATS is a function that performs statistical analyses to compare
    U-Pb age probability density functions (PDFs).

     Description
    
       'zSTATS' is a function that performs 4 different statistical comparison
       of U-Pb age populations. It takes the PDF names and sizes, and returns
       similarity, area mismatch, likeliness and Kolmogorov-Smirnov test.

       For a detailed explanation of the equations, please refer to the work of
       Lavarini et al. [2018a].
    
     Input arguments
    
        Ami = PDF 1
        Lav = PDF 2 
          
     Output arguments 
    
       It returns printed values of similarity, area mismatch, likeliness 
       and Kolmogorov-Smirnov test results.
    
     Reference:
         
         Lavarini et al. [2018]. Does pebble abrasion influence detrital age 
         population statistics? A numerical investigation of natural data sets. 
         Journal of Geophysical Research - Earth Surface.
         
         Amidon et al. [2005a]. Construction of detrital mineral populations: 
             Insights from mixing of U-Pb zircon ages in Himalayan rivers.
             Basin Research.
    
     Authors: C. Lavarini (lavarini.c@gmail.com), C. A. da Costa Filho & M. Attal
     The University of Edinburgh, UK
     December of 2017
    
    """
    # Similarity
    S1 = np.sqrt((Ami)*(Lav)).sum()
    # Area mismatch
    M1 = np.abs((Ami)-(Lav)).sum()/2
    # Likeliness
    L1 = 1-M1  
    
    # Display the stats:
    print('S1 =', S1)    
    print('-------')
    print('M1 =', M1)    
    print('-------')
    print('L1 =', L1)
    
    # Kolmogorov-Smirnov



x1 = [Lav.as_matrix().flatten(), Ami.as_matrix().flatten()]

# Q-Q plot

def qqPLOT(x1):

    """
   
    Q-Q plot  ('qqPLOT') is a function that generates cumulative distribution 
    functions (CDFs) of U-Pb ages.
    
     Description
    
       'qqPLOT' is a function that reads PDFs of U-Pb ages and then 
       creates a regression analysis in cumulative distribution 
       functions (CDFs)for a pair of U-Pb age list.
       
       For a detailed explanation of the equations, please refer to the work of
       Lavarini et al. [2018a].

     Input arguments
    
       x1 = A list with a pair of PDFs to be compared.
                  
     Output arguments 
    
       It returns a Q-Q plot of the two samples.
             
     Example
    
    x1 = [Lav.as_matrix().flatten(), Ami.as_matrix().flatten()]

     Reference:
         
         Lavarini et al. [2018]. Does pebble abrasion influence detrital age 
         population statistics? A numerical investigation of natural data sets. 
         Journal of Geophysical Research - Earth Surface.
    
     Authors: C. Lavarini (lavarini.c@gmail.com), C. A. da Costa Filho & M. Attal
     The University of Edinburgh, UK
     December of 2017
    
    """
    # Create CDFs
    dx = Lav.index[1] - Lav.index[0] 
    x = cumtrapz(Lav.as_matrix().flatten(), dx=dx)
    y = cumtrapz(Ami.as_matrix().flatten(), dx=dx)
    
    axmin, axmax = min([min(x), min(y)]), max([max(x), max(y)])
    axmin,axmax = axmin-(axmax-axmin)*0.05, axmax+(axmax-axmin)*0.05

    X = sm.add_constant(x)
    reg0 = sm.OLS(y, X).fit()
   
    f, ax = plt.subplots()
    ax.plot(x, x, 'k')
    
    ax.scatter(x, y,  s=50, alpha=0.6, edgecolor='r', facecolors='none', label=r'$R^2 = %.3f$' % reg0.rsquared)
    
    ax.set_xlabel('Artificial CDF [K]')
    ax.set_ylabel('Predicted CDFs')
    ax.set_title('Q-Q plot')
    #sns.despine(ax=ax)
    ax.set_xlim([axmin, axmax])
    ax.set_ylim([axmin, axmax])
    plt.legend(loc=4)
    handles, labels = ax.get_legend_handles_labels()
    leg = ax.legend(handles, labels, frameon=False, fancybox=True, fontsize=18)
    leg.get_frame().set_edgecolor('k')
    f.savefig(figpath+'crossplotcdf_G.png', bbox_inches='tight')
    f.savefig(figpath+'crossplotcdf_G.pdf', bbox_inches='tight')

x1 = [Lav.as_matrix().flatten(), Ami.as_matrix().flatten()]

# PDF crossplot
def pdfCROSS(x1):
    
    """
   
    PDF crossplot  ('pdfCROSS') is a function that generates probability 
    density function (PDFs) plots of U-Pb ages.
    
     Description
    
       'pdfCROSS' is a function that reads PDFs of U-Pb ages and then 
       creates a regression analysis in probability density plots (PDPs) for a 
       pair of U-Pb age list.
       
       For a detailed explanation of the equations, please refer to the work of
       Lavarini et al. [2018a].

     Input arguments
    
       x1 = A list with a pair of PDFs to be compared.
                  
     Output arguments 
    
       It returns the PDF crossplot of the two samples.
             
     Example
    
    x1 = [Lav.as_matrix().flatten(), Ami.as_matrix().flatten()]

     Reference:
         
         Lavarini et al. [2018]. Does pebble abrasion influence detrital age 
         population statistics? A numerical investigation of natural data sets. 
         Journal of Geophysical Research - Earth Surface.
    
     Authors: C. Lavarini (lavarini.c@gmail.com), C. A. da Costa Filho & M. Attal
     The University of Edinburgh, UK
     December of 2017
    
    """
    
    X = sm.add_constant(x1[0])
    reg1 = sm.OLS(x1[1], X).fit()
   
    sns.set(style='ticks', font_scale=1.5)
    
    # Pretty
    #cm = plt.cm.get_cmap('RdYlBu')
    f, ax = plt.subplots()
    ax.plot(x1[0], x1[0], 'k')
    ax.scatter(x1[0], x1[1],  s=50, alpha=0.6, edgecolor='b', \
               facecolors='none', label=r'$R^2 = %.3f$' % reg1.rsquared)

    ax.set_title('PDF cross-plots')
    ax.set_xlabel('PDF')
    ax.set_ylabel('PDF')
    plt.ticklabel_format(style='sci', scilimits=(-2,2))   

    plt.legend(loc=4)
    handles, labels = ax.get_legend_handles_labels()
    leg = ax.legend(handles, labels, frameon=False, fancybox=True, fontsize=18)
    leg.get_frame().set_edgecolor('k')
    print('PDF crossplot generated!')    