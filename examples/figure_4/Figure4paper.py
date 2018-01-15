#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
 
Figure 4 (Figure4paper.py) is a series of code used to generate the figure 4
in the work of Lavarini et al. [2018a].
 Reference:
     
     Lavarini et al. [2018]. Does pebble abrasion influence detrital age 
     population statistics? A numerical investigation of natural data sets. 
     Journal of Geophysical Research - Earth Surface.
     
 Authors: C. Lavarini (lavarini.c@gmail.com), C. A. da Costa Filho & M. Attal
 The University of Edinburgh, UK
 December of 2017
      
"""

# Import modules
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from numpy import pi, exp, sqrt
from scipy import stats
from scipy.optimize import minimize
from scipy.integrate import cumtrapz
import statsmodels.api as sm
blu, grn, red, prp, yel, cyn = sns.color_palette()

# Set the style
sns.set(style='ticks', font_scale=1.5, 
        rc={'lines.linewidth': 2.5,
            'figure.figsize': (10, 8),
            'text.usetex': False,
        })
# Set path to save figures:
figpath = os.path.join('K_AGES/', 'marsyandi_')

# Function to generate probability density functions (PDFs):
def generate_pdf(filename, sheet, smooth=80, plot=True):
    df = pd.read_excel(filename, sheetname = sheet, skiprows=(0,))

    # Create single-column with the true age
    K = df[['Age (Ma)', 'Error (Ma)']]
    K.columns = ['m', 's']
    K.name = 'Sample '+sheet

    K['m'] = df['Age (Ma)'] 
    K['s']  = df['Error (Ma)'] 

    # Show Histogram
    if plot:
        f, ax = plt.subplots()
        plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        K['m'].hist(bins=50, ax=ax, color='0.6')
        ax.set_xlim([0, 4000])
        ax.set_xlabel('Synthetic U-Pb age [Ma]')
        ax.set_ylabel('Frequency')
        f.savefig(figpath+'histo_'+sheet+'.pdf', bbox_inches='tight')
        f.savefig(figpath+'histo_'+sheet+'.png', bbox_inches='tight')
    sample_size = len(K['m'])
    
    # Calculate pdf
    x = np.linspace(0, 4000, 4001) # Ma
    pdf = np.zeros(len(x))
    for i, row in K.iterrows():
        # Create gaussian for each sample
        mu_i, sigma_i = row['m'], row['s']
        pdf += exp( -(x - mu_i)**2/(2*sigma_i**2) )/(2*sigma_i*sqrt(2*pi))

    if plot:
        f, ax = plt.subplots()
        plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        ax.fill_between(x, 0, pdf, color='k')
        ax.set_xlabel('Synthetic U-Pb age [Ma]')
        ax.set_ylabel('Relative probability')
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
        f.savefig(figpath+'pdf_smooth_'+sheet+'.pdf', bbox_inches='tight')
        f.savefig(figpath+'pdf_smooth_'+sheet+'.png', bbox_inches='tight')
    return pdf, sample_size

#######################
## Generate all PDFs ##
#######################

plt.ion()

# Sheets with single source unit data
sheets = ['A', 'C', 'F', 'H', 'K']

# Filename containing the sheets
filename = 'MARSYANDI_SYNTHETIC_1.xlsx'

# Sheet names of single source units
objects = ('A', 'C', 'K', 'H', 'F')
objects2 = ('A', 'C', 'H', 'F')

# Create a dictionary to store the PDFs and their sample sizes
PDF = {}
SIZES = {}

# Loop over all the sheets (samples)
for sheet in sheets:
    PDF[sheet], SIZES[sheet] = generate_pdf(filename, sheet, smooth=80, plot=False)

# Testing the influence of controlling factors in zircon mixing proportion

# No factors occurring. 
# Values of mixing proportions from pABRASIONmodel where 0 %/km was used for all sources (i.e., no abrasion case)
phi_no = [0.271, 0.175, 0.235, 0.319]  
# Create PDF
PDF['Ami'] = phi_no[0]*PDF['A'] + phi_no[1]*PDF['C'] + \
         phi_no[2]*PDF['F'] + phi_no[3]*PDF['H'] 
 
# Erosion.
# Values of mixing proportions from pABRASIONmodel where TTS erodes at 5.1 mm/yr and the rest at 1.0 mm/yr.
phi_ami2 = [0.654, 0.083, 0.111, 0.152]
# Create PDF
PDF['Ami2'] = phi_ami2[0]*PDF['A'] + phi_ami2[1]*PDF['C'] + \
         phi_ami2[2]*PDF['F'] + phi_ami2[3]*PDF['H']
 
# No factors occurring. 
# Values of mixing proportions from pABRASIONmodel where 0 %/km was used for all sources (i.e., no abrasion case)
# Create PDF
PDF['Abr'] = phi_no[0]*PDF['A'] + phi_no[1]*PDF['C'] + \
         phi_no[2]*PDF['F'] + phi_no[3]*PDF['H'] 

# Zircon fertility 
# Values of mixing proportions from pABRASIONmodel where 8.1 grains/g for TTS and 0 grains/g for rest
# Create PDF
PDF['Art'] = 1*PDF['A'] + 0*PDF['C'] + 0*PDF['F'] + 0*PDF['H'] 

# No factors occurring. 
PDF['Art3'] = phi_no[0]*PDF['A'] + phi_no[1]*PDF['C'] + \
         phi_no[2]*PDF['F'] + phi_no[3]*PDF['H'] # No variables

# Hillslope gravel supply
# Values of mixing proportions from pABRASIONmodel where 50% of sediments from TTS is sand and 15% of the remaining 
# watershed supplies sand.
# Create PDF
PDF['Art31'] = 0.565*PDF['A'] + 0.104*PDF['C'] + 0.14*PDF['F'] + 0.191*PDF['H'] 

# Pebble abrasion rate
# Values of mixing proportions from pABRASIONmodel where 31 %/km was used for TTS
# and 0.15  %/km for all other sources (i.e., no abrasion case)
phi_abr = [0.525, 0.129, 0.158, 0.188] 
# Create PDF
PDF['Art41'] = phi_abr[0]*PDF['A'] + phi_abr[1]*PDF['C'] + \
         phi_abr[2]*PDF['F'] + phi_abr[3]*PDF['H'] 
 
# No factors occurring. 
PDF['Art4'] = phi_no[0]*PDF['A'] + phi_no[1]*PDF['C'] + \
         phi_no[2]*PDF['F'] + phi_no[3]*PDF['H'] 

# Erosion vs no factor
# Similarity:
S1 = np.sqrt(PDF['Ami']*PDF['Ami2']).sum()
# Mismatch:
M1 = np.abs(PDF['Ami']-PDF['Ami2']).sum()/2
# Likeliness:
L1 = 1-M1

# Fertility vs no factor
# Similarity:
S2 = np.sqrt(PDF['Abr']*PDF['Art']).sum()
M2 = np.abs(PDF['Abr']-PDF['Art']).sum()/2
L2 = 1-M2

# Gravel supply vs no factor
# Similarity:
S3 = np.sqrt(PDF['Art3']*PDF['Art31']).sum()
# Mismatch:
M3 = np.abs(PDF['Art3']-PDF['Art31']).sum()/2
# Likeliness:
L3 = 1-M3

# Abrasion rate vs no factor
# Similarity:
S4 = np.sqrt(PDF['Art4']*PDF['Art41']).sum()
# Mismatch:
M4 = np.abs(PDF['Art4']-PDF['Art41']).sum()/2
# Likeliness:
L4 = 1-M4

#  Display calculations:
print('S1 =', S1)
print('S2 =', S2)
print('S3 =', S3)
print('S4 =', S4)
print('-------')
print('M1 =', M1)
print('M2 =', M2)
print('M3 =', M3)
print('M4 =', M4)
print('-------')
print('L1 =', L1)
print('L2 =', L2)
print('L3 =', L3)
print('L4 =', L4)


# Distance stored
x = PDF['Ami'].index

# Create distribution for Kolmogorv-Smirnov (K-S) test
DIST = {}

DIST['Ami'] = stats.rv_discrete(name='Ami',
                values=(PDF['Ami'].index, PDF['Ami'].as_matrix().flatten()))
DIST['Ami2'] = stats.rv_discrete(name='Ami2',
                values=(PDF['Ami2'].index, PDF['Ami2'].as_matrix().flatten()))
DIST['Abr'] = stats.rv_discrete(name='Abr',
                values=(PDF['Abr'].index, PDF['Abr'].as_matrix().flatten()))
DIST['Art'] = stats.rv_discrete(name='Art',
                values=(PDF['Art'].index, PDF['Art'].as_matrix().flatten()))
DIST['Art3'] = stats.rv_discrete(name='Art3',
                values=(PDF['Art3'].index, PDF['Art3'].as_matrix().flatten()))
DIST['Art31'] = stats.rv_discrete(name='Art31',
                values=(PDF['Art3'].index, PDF['Art31'].as_matrix().flatten()))
DIST['Art4'] = stats.rv_discrete(name='Art4',
                values=(PDF['Art4'].index, PDF['Art4'].as_matrix().flatten()))
DIST['Art41'] = stats.rv_discrete(name='Art41',
                values=(PDF['Art4'].index, PDF['Art41'].as_matrix().flatten()))

# Iterate over all single source data
for region in objects:
    DIST[region] = stats.rv_discrete(name=region,
                                values=(x, PDF[region].as_matrix().flatten()))

# Calculate KS test on two sampled distributions
SIZES['Art'] = SIZES['A'] + SIZES['F'] + SIZES['C'] + SIZES['H']
D1, p1 = 0, 0
D2, p2 = 0, 0
D3, p3 = 0, 0
D4, p4 = 0, 0
D5, p5 = 0, 0

# Storing the test results
for i in range(1000):
    data1, data2 = DIST['Ami'].rvs(size=SIZES['Art']), DIST['Ami2'].rvs(size=SIZES['Art'])
    data3, data4 = DIST['Abr'].rvs(size=SIZES['Art']), DIST['Art'].rvs(size=SIZES['Art'])
    data5, data6 = DIST['Art3'].rvs(size=SIZES['Art']), DIST['Art31'].rvs(size=SIZES['Art'])
    data7, data8 = DIST['Art4'].rvs(size=SIZES['Art']), DIST['Art41'].rvs(size=SIZES['Art'])
  
    kstest1 = stats.ks_2samp(data1, data2)
    kstest2 = stats.ks_2samp(data3, data4)
    kstest3 = stats.ks_2samp(data5, data6)
    kstest4 = stats.ks_2samp(data7, data8)
   
    D1 += kstest1.statistic
    p1 += kstest1.pvalue
    D2 += kstest2.statistic
    p2 += kstest2.pvalue
    D3 += kstest3.statistic
    p3 += kstest3.pvalue
    D4 += kstest4.statistic
    p4 += kstest4.pvalue
    
D1 /= 1000
p1 /= 1000
D2 /= 1000
p2 /= 1000
D3 /= 1000
p3 /= 1000
D4 /= 1000
p4 /= 1000

# Display the results:
print('-------')
print('KS 1', kstest1)
print('KS 2', kstest2)
print('KS 3', kstest3)
print('KS 4', kstest4)

print('-------')
print('KS 1a', D1, p1)
print('KS 2a', D2, p2)
print('KS 3a', D3, p3)
print('KS 4a', D4, p4)
print('-------')

# Plot the probability density functions (PDFs) generated in the previous scenarios:
f, ax = plt.subplots()
plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
# PDF of a no factor case.
ax.plot(x, PDF['Ami'], label='Observed', color = '0.6')
ax.fill_between(x, 0, PDF['Ami'].as_matrix().flatten(), color = '0.6')
# PDF of erosion scenario:
ax.plot(x, PDF['Ami2'], 'b',
           label='1b')
# PDF of fertility scenario:
ax.plot(x, PDF['Art'], 'g',
           label='2b')
# PDF of gravel scenario:
ax.plot(x, PDF['Art31'], 'y',
           label='3b')
# PDF of abrasion scenario:
ax.plot(x, PDF['Art41'], 'c',
           label='4b')
# Legend:
handles, labels = ax.get_legend_handles_labels()
leg = ax.legend(handles, labels, frameon=True, fancybox=True, fontsize=18)
leg.get_frame().set_edgecolor('k')
ax.set_xlabel('Age [Ma]')
ax.set_ylabel('Relative probability')
ax.set_title('PDPs')
# Save plot:
f.savefig(figpath+'pdf_mix_K_fig_4.pdf', bbox_inches='tight')
f.savefig(figpath+'pdf_mix_K_fig_4.png', bbox_inches='tight')

# Create CDFs for Q-Q plots
# Measuring dx:
dx = PDF['K'].index[1] - PDF['K'].index[0] 
# CDF of no-abrasion case:
xa = cumtrapz(PDF['Ami'].as_matrix().flatten(), dx=dx)
# CDF of erosion scenario:
xb = cumtrapz(PDF['Ami2'].as_matrix().flatten(), dx=dx)
# CDF of no-abrasion case:
za = cumtrapz(PDF['Abr'].as_matrix().flatten(), dx=dx)
# CDF of fertility scenario:
zb = cumtrapz(PDF['Art'].as_matrix().flatten(), dx=dx)
# CDF of no-abrasion case:
ya = cumtrapz(PDF['Art3'].as_matrix().flatten(), dx=dx)
# CDF of gravel scenario:
yb = cumtrapz(PDF['Art31'].as_matrix().flatten(), dx=dx)
# CDF of no-abrasion case:
ka = cumtrapz(PDF['Art4'].as_matrix().flatten(), dx=dx)
# CDF of abrasion scenario:
kb = cumtrapz(PDF['Art41'].as_matrix().flatten(), dx=dx)

# Max and min of values
axmin, axmax = min([min(xa), min(xb)]), max([max(xa), max(xb)])
axmin,axmax = axmin-(axmax-axmin)*0.05, axmax+(axmax-axmin)*0.05

# Set the Y axis:
Xa = sm.add_constant(xb)
Xb = sm.add_constant(zb)
Xc = sm.add_constant(yb)
Xd = sm.add_constant(kb)

# Set the X axis for regression:
reg0 = sm.OLS(xa, Xa).fit()
reg1 = sm.OLS(za, Xb).fit()
reg2 = sm.OLS(ya, Xc).fit()
reg3 = sm.OLS(ka, Xd).fit()

# Plot the CDFs:
f, ax = plt.subplots()
ax.plot(xa, xb, 'k')

# CDF of erosion scenario:
ax.scatter(xa, xb,  s=50, alpha=0.6, edgecolor='b', facecolors='none', label=r'$R^2 = %.3f$' % reg0.rsquared)
# CDF of fertility scenario:
ax.scatter(za, zb, s=50, alpha=0.6, edgecolor='g', facecolors='none', label=r'$R^2 = %.3f$' % reg1.rsquared)
# CDF of gravel scenario:
ax.scatter(ya, yb, s=50, alpha=0.6, edgecolor='y', facecolors='none',label=r'$R^2 = %.3f$' % reg2.rsquared)
# CDF of abrasion scenario:
ax.scatter(ka, kb, s=50, alpha=0.6, edgecolor='c', facecolors='none',label=r'$R^2 = %.3f$' % reg3.rsquared)
# Set the legends:
ax.set_xlabel('CDFs [Tested factors]')
ax.set_ylabel('CDFs [Target factors]')
ax.set_title('Q-Q plots')
ax.set_xlim([axmin, axmax])
ax.set_ylim([axmin, axmax])
plt.legend(loc=4)
leg.get_frame().set_edgecolor('k')
# Save the plot:
f.savefig(figpath+'crossplotcdf_G.png', bbox_inches='tight')
f.savefig(figpath+'crossplotcdf_G.pdf', bbox_inches='tight')

# Making PDF cross-plots
# List with a PDF for each scenario:
# PDFs made with no factor
x1 = [PDF['Ami'].as_matrix().flatten(), PDF['Abr'].as_matrix().flatten(), PDF['Art3'].as_matrix().flatten(), PDF['Art4'].as_matrix().flatten()]
# PDFs made with testing factors:
x2 = [PDF['Ami2'].as_matrix().flatten(), PDF['Art'].as_matrix().flatten(),PDF['Art31'].as_matrix().flatten(), PDF['Art41'].as_matrix().flatten()]
# Set the Y axis:
X1 = sm.add_constant(x2[0])
X2 = sm.add_constant(x2[1])
X3 = sm.add_constant(x2[2])
X4 = sm.add_constant(x2[3])
# Set the X axis for regression:
reg1 = sm.OLS(x1[0], X1).fit()
reg2 = sm.OLS(x1[1], X2).fit()
reg3 = sm.OLS(x1[2], X3).fit()
reg4 = sm.OLS(x1[3], X4).fit()

# Set the plot style:
sns.set(style='ticks', font_scale=1.5)

# Plot the PDF crossplots:
f, ax = plt.subplots()
# PDF of a no-abrasion case:
ax.plot(x1[0], x2[0])
# PDF of erosion scenario:
ax.scatter(x1[0], x2[0],  s=50, alpha=0.6, edgecolor='b', facecolors='none', label=r'$R^2 = %.3f$' % reg1.rsquared)
# PDF of fertility scenario:
ax.scatter(x1[1], x2[1], s=50, alpha=0.6, edgecolor='g', facecolors='none', label=r'$R^2 = %.3f$' % reg2.rsquared)
# PDF of gravel scenario:
ax.scatter(x1[2], x2[2], s=50, alpha=0.6, edgecolor='y', facecolors='none',label=r'$R^2 = %.3f$' % reg3.rsquared)
# PDF of abrasion scenario:
ax.scatter(x1[3], x2[3], s=50, alpha=0.6, edgecolor='c', facecolors='none',label=r'$R^2 = %.3f$' % reg4.rsquared)
# Set the legends:
ax.set_title('PDF cross-plots')
ax.set_xlabel('Tested PDFs')
ax.set_ylabel('Target PDFs')
plt.ticklabel_format(style='sci', scilimits=(-2,2))   
plt.legend(loc=4)
leg.get_frame().set_edgecolor('k')
# Save the plot:
f.savefig(figpath+'crossplotpdf_K2.png', bbox_inches='tight')
f.savefig(figpath+'crossplotpdf_K2.pdf', bbox_inches='tight')


