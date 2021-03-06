#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 
Figure 3 (Figure3paper.py) is a series of code used to generate the figure 3
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

# Function to find the best-fit zircon mixing proportion
def cost_function(phi):
     f = phi[0]*PDF['A'] + phi[1]*PDF['C'] + \
         phi[2]*PDF['F'] + phi[3]*PDF['H']
     x = f.index
     cost = np.trapz(abs(f-PDF['K']).as_matrix().flatten(),
                      dx=x[2]-x[1])/2
     return cost*100
    
# Function to find the best-fit zircon mixing proportion
def cost_function2(phi):
     f = phi[0]*0.31*PDF['A'] + phi[1]*0.12*PDF['C'] + \
         phi[2]*0.23*PDF['F'] + phi[3]*0.33*PDF['H']
     x = f.index
     cost = np.trapz(abs(f-PDF['K']).as_matrix().flatten(),
                      dx=x[2]-x[1])/2
     return cost*100
    
# Finding the best-fit zircon mixing proportion
initial_phi = [0, 1, 0, 0]
bounds = tuple([(0,1)]*len(initial_phi))
constraints = ({'type':'eq', 'fun': lambda x: sum(x)-1},)
res = minimize(cost_function, initial_phi, method='SLSQP', bounds=bounds, constraints=constraints, tol=1e-10)
lowest_cost = res.fun
best_phi = res.x
for i in range(30):
    initial_phi = np.random.uniform(0, 1, size=(len(initial_phi),))
    initial_phi /= sum(initial_phi)
    bounds = tuple([(0,1)]*len(initial_phi))
    res = minimize(cost_function, initial_phi, method='SLSQP', bounds=bounds, constraints=constraints, tol=1e-10)
    if res.fun < lowest_cost:
        lowest_cost = res.fun
        best_phi = res.x
        print()
        print('New best cost: ', lowest_cost)
        print('New best phi:  ', best_phi)
        
        
# Scenario A1
# Testing the influence of abrasion in zircon mixing proportion
# Values of mixing proportions from pABRASIONmodel where 0.4%/km was used for all sources
phi_ami = [0.331, 0.194, 0.229, 0.246] 
# Generate PDFs
PDF['Ami'] = phi_ami[0]*PDF['A'] + phi_ami[1]*PDF['C'] + \
         phi_ami[2]*PDF['F'] + phi_ami[3]*PDF['H']
x = PDF['Ami'].index
# Calculate area mismatch between the above PDF and a PDF with no-abrasion
newcost = np.trapz(abs(PDF['Ami']-PDF['K']).as_matrix().flatten(),
                      dx=x[2]-x[1])/2
newcost = newcost*100
# Print the mismatch and zircon mixing proportion:
print('-------')
print('Mismatch:', newcost)
print('Zircon mixing proportion:', phi_ami)


# Scenario A2
# Testing the influence of abrasion in zircon mixing proportion
# Values of mixing proportions from pABRASIONmodel where 31.0 %/km was used for the TTS and 0.15 for all rest
phi_ami2 = [0.525, 0.129, 0.158, 0.188] 
# Generate PDFs
PDF['Ami2'] = phi_ami2[0]*PDF['A'] + phi_ami2[1]*PDF['C'] + \
         phi_ami2[2]*PDF['F'] + phi_ami2[3]*PDF['H']
x = PDF['Ami2'].index
# Calculate area mismatch between the above PDF and a PDF with no-abrasion
newcost2 = np.trapz(abs(PDF['Ami2']-PDF['K']).as_matrix().flatten(),
                      dx=x[2]-x[1])/2
newcost2 = newcost2*100
# Print the mismatch and zircon mixing proportion:
print('-------')
print('Mismatch:', newcost2)
print('Zircon mixing proportion:', phi_ami2)


# Scenario A3
# Testing the influence of abrasion in zircon mixing proportion
# Values of mixing proportions from pABRASIONmodel where 31.0 %/km was used for the LH and 0.15 for all rest
phi_abr = [0.192, 0.115, 0.140, 0.553] 
# Generate PDFs
PDF['Abr'] = phi_abr[0]*PDF['A'] + phi_abr[1]*PDF['C'] + \
         phi_abr[2]*PDF['F'] + phi_abr[3]*PDF['H']
g = PDF['Abr'].index
# Calculate area mismatch between the above PDF and a PDF with no-abrasion
newcost3 = np.trapz(abs(PDF['Abr']-PDF['K']).as_matrix().flatten(),
                      dx=x[2]-x[1])/2
newcost3 = newcost3*100
# Print the mismatch and zircon mixing proportion:
print('-------')
print('Mismatch:', newcost3)
print('Zircon mixing proportion:', phi_abr)

# Scenario A4
# Testing the influence of abrasion in zircon mixing proportion
# Values of mixing proportions from pABRASIONmodel where the values are the same found by Attal and Lave [2006].
PDF['Art'] = 0.318*PDF['A'] + 0.112*PDF['C'] + 0.207*PDF['F'] + 0.363*PDF['H'] 


# Statistical comparison of the PDFs generated in the scenarios above (A1 to A4):

# Scenario A1
# Similarity
S1 = np.sqrt(PDF['Ami']*PDF['K']).sum()
# Mismatch
M1 = np.abs(PDF['Ami']-PDF['K']).sum()/2
# Likeness
L1 = 1-M1
# Display results
print('--------')
print('S1 =', S1)
print('M1 =', M1)
print('L1 =', L1)
print('--------')

# Scenario A2
# Similarity
S2 = np.sqrt(PDF['Ami2']*PDF['K']).sum()
# Mismatch
M2 = np.abs(PDF['Ami2']-PDF['K']).sum()/2
# Likeness
L2 = 1-M2
# Display results
print('--------')
print('S2 =', S2)
print('M2 =', M2)
print('L2 =', L2)
print('--------')

# Scenario A3
# Similarity
S3 = np.sqrt(PDF['Abr']*PDF['K']).sum()
# Mismatch
M3 = np.abs(PDF['Abr']-PDF['K']).sum()/2
# Likeness
L3 = 1-M3
# Display results
print('--------')
print('S3 =', S3)
print('M3 =', M3)
print('L3 =', L3)
print('--------')

# Scenario A4
# Similarity
S4 = np.sqrt(PDF['Art']*PDF['K']).sum()
# Mismatch
M4 = np.abs(PDF['Art']-PDF['K']).sum()/2
# Likeness
L4 = 1-M4
# Display results
print('--------')
print('S4 =', S4)
print('M4 =', M4)
print('L4 =', L4)
print('--------')


# Create distribution for Kolmogorv-Smirnov (K-S) test
DIST = {}
DIST['Art'] = stats.rv_discrete(name='Art',
                values=(x, PDF['Art'].as_matrix().flatten()))
DIST['K'] = stats.rv_discrete(name='K',
                values=(x, PDF['K'].as_matrix().flatten()))
DIST['Ami'] = stats.rv_discrete(name='Ami',
                values=(x, PDF['Ami'].as_matrix().flatten()))
DIST['Ami2'] = stats.rv_discrete(name='Ami2',
                values=(x, PDF['Ami2'].as_matrix().flatten()))
DIST['Abr'] = stats.rv_discrete(name='Abr',
                values=(x, PDF['Abr'].as_matrix().flatten()))

# Iterate over all single source data
for region in objects:
    DIST[region] = stats.rv_discrete(name=region,
                                values=(x, PDF[region].as_matrix().flatten()))

# Calculate KS test on two sampled distributions
SIZES['Art'] = SIZES['A'] + SIZES['F'] + SIZES['C'] + SIZES['H']
# Distance between two CDFs (D) and its probability (p) of rejecting the null hypothesis:
D1, p1 = 0, 0
D2, p2 = 0, 0
D3, p3 = 0, 0
D4, p4 = 0, 0

# Storing the test results
for i in range(1000):
    data1, data2 = DIST['Art'].rvs(size=SIZES['Art']), DIST['K'].rvs(size=SIZES['K'])
    data3, data4 = DIST['Ami'].rvs(size=SIZES['Art']), DIST['K'].rvs(size=SIZES['K'])
    data5, data6 = DIST['Ami2'].rvs(size=SIZES['Art']), DIST['K'].rvs(size=SIZES['K'])
    data7, data8 = DIST['Abr'].rvs(size=SIZES['Art']), DIST['K'].rvs(size=SIZES['K'])
    
    kstest1 = stats.ks_2samp(data3, data4)
    kstest2 = stats.ks_2samp(data5, data6)
    kstest3 = stats.ks_2samp(data7, data8)
    kstest4 = stats.ks_2samp(data1, data2)
    
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
print('KS 1: ', kstest1)
print('KS 2: ', kstest2)
print('KS 3: ', kstest3)
print('KS 4: ', kstest4)
print('-------')
print('KS 1: ', D1, p1)
print('KS 2: ', D2, p2)
print('KS 3: ', D3, p3)
print('KS 4: ', D4, p4)
print('-------')


# Plot the probability density functions (PDFs) generated in the previous scenarios (A1-A4):
f, ax = plt.subplots()
plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
# PDF of a no-abrasion case.
ax.plot(x, PDF['K'], label='Observed', color = '0.6')
ax.fill_between(x, 0, PDF['K'].as_matrix().flatten(), color = '0.6')
# PDF of scenario A1:
ax.plot(x, PDF['Ami'], 'b',
           label='1')
# PDF of scenario A2:
ax.plot(x, PDF['Ami2'], 'y',
           label='2')
# PDF of scenario A3:
ax.plot(x, PDF['Abr'], 'g',
           label='3')
# PDF of scenario A4:
ax.plot(x, PDF['Art'], 'r',
           label='4')
handles, labels = ax.get_legend_handles_labels()
# Legend:
leg = ax.legend(handles, labels, frameon=True, fancybox=True, fontsize=18)
leg.get_frame().set_edgecolor('k')
ax.set_xlabel('Age [Ma]')
ax.set_ylabel('Relative probability')
ax.set_title('PDPs (X)')
# Save plot:
f.savefig(figpath+'pdf_mix_K_E.pdf', bbox_inches='tight')
f.savefig(figpath+'pdf_mix_K_E.png', bbox_inches='tight')

# Create CDFs for Q-Q plot
# Measuring dx:
dx = PDF['K'].index[1] - PDF['K'].index[0] 
# CDF of no-abrasion case:
x = cumtrapz(PDF['K'].as_matrix().flatten(), dx=dx)
# CDF of scenario A1:
y = cumtrapz(PDF['Art'].as_matrix().flatten(), dx=dx)
# CDF of scenario A2:
xh = cumtrapz(PDF['Ami'].as_matrix().flatten(), dx=dx)
# CDF of scenario A3:
xj = cumtrapz(PDF['Ami2'].as_matrix().flatten(), dx=dx)
# CDF of scenario A4:
jj = cumtrapz(PDF['Abr'].as_matrix().flatten(), dx=dx)
# Max and min of values
axmin, axmax = min([min(x), min(y)]), max([max(x), max(y)])
axmin,axmax = axmin-(axmax-axmin)*0.05, axmax+(axmax-axmin)*0.05

# Set the Y axis:
X = sm.add_constant(x)
# Set the X axis for regression:
reg0 = sm.OLS(y, X).fit()
reg1 = sm.OLS(xh, X).fit()
reg2 = sm.OLS(xj, X).fit()
reg3 = sm.OLS(jj, X).fit()

# Plot the CDFs:
f, ax = plt.subplots()
ax.plot(x, x, 'k')

# CDF of scenario A1:
ax.scatter(x, y,  s=50, alpha=0.6, edgecolor='r', facecolors='none', label=r'$R^2 = %.3f$' % reg0.rsquared)
# CDF of scenario A2:
ax.scatter(x, xh, s=50, alpha=0.6, edgecolor='b', facecolors='none', label=r'$R^2 = %.3f$' % reg1.rsquared)
# CDF of scenario A3:
ax.scatter(x, xj, s=50, alpha=0.6, edgecolor='y', facecolors='none',label=r'$R^2 = %.3f$' % reg2.rsquared)
# CDF of scenario A4:
ax.scatter(x, jj, s=50, alpha=0.6, edgecolor='g', facecolors='none', label=r'$R^2 = %.3f$' % reg3.rsquared)
# Set the legends:
ax.set_xlabel('Artificial CDF [K]')
ax.set_ylabel('Predicted CDFs')
ax.set_title('Q-Q plot')
ax.set_xlim([axmin, axmax])
ax.set_ylim([axmin, axmax])
plt.legend(loc=4)
leg.get_frame().set_edgecolor('k')
# Save the plot:
f.savefig(figpath+'crossplotcdf_K.png', bbox_inches='tight')
f.savefig(figpath+'crossplotcdf_K.pdf', bbox_inches='tight')

# PDF crossplots:
# List with a PDF for each scenario (A1 to A4):
x1 = [PDF['K'].as_matrix().flatten(), PDF['Ami'].as_matrix().flatten(), PDF['Ami2'].as_matrix().flatten(), PDF['Abr'].as_matrix().flatten(), PDF['Art'].as_matrix().flatten()]

# Set the Y axis:
X = sm.add_constant(x1[0])
# Set the X axis for regression:
reg1 = sm.OLS(x1[1], X).fit()
reg2 = sm.OLS(x1[2], X).fit()
reg3 = sm.OLS(x1[3], X).fit()
reg4 = sm.OLS(x1[4], X).fit()

# Set the plot style:
sns.set(style='ticks', font_scale=1.5)

# Plot the PDF crossplots:
f, ax = plt.subplots()

# PDF of a no-abrasion case:
ax.plot(x1[0], x1[0], 'k')
# PDF of scenario A1:
ax.scatter(x1[0], x1[1],  s=50, alpha=0.6, edgecolor='b', facecolors='none', label=r'$R^2 = %.3f$' % reg1.rsquared)
# PDF of scenario A2:
ax.scatter(x1[0], x1[2], s=50, alpha=0.6, edgecolor='y', facecolors='none', label=r'$R^2 = %.3f$' % reg2.rsquared)
# PDF of scenario A3:
ax.scatter(x1[0], x1[3], s=50, alpha=0.6, edgecolor='g', facecolors='none',label=r'$R^2 = %.3f$' % reg3.rsquared)
# PDF of scenario A4:
ax.scatter(x1[0], x1[4], s=50, alpha=0.6, edgecolor='r', facecolors='none', label=r'$R^2 = %.3f$' % reg4.rsquared)
# Set the legends:
ax.set_title('PDF cross-plots')
ax.set_xlabel('PDF')
ax.set_ylabel('PDF')
plt.ticklabel_format(style='sci', scilimits=(-2,2))   
plt.legend(loc=4)
leg.get_frame().set_edgecolor('k')
# Save the plot:
f.savefig(figpath+'crossplotpdf_K.png', bbox_inches='tight')
f.savefig(figpath+'crossplotpdf_K.pdf', bbox_inches='tight')
