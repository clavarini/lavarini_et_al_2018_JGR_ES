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
print('Lowest mismatch Amidon paper:', newcost2)
print('Best phi Amidon paper:', phi_ami2)

phi_abr = [0.192, 0.115, 0.140, 0.553] # 31 for the last. 0.15 for all rest

PDF['Abr'] = phi_abr[0]*PDF['A'] + phi_abr[1]*PDF['C'] + \
         phi_abr[2]*PDF['F'] + phi_abr[3]*PDF['H']
g = PDF['Abr'].index
newcost3 = np.trapz(abs(PDF['Abr']-PDF['K']).as_matrix().flatten(),
                      dx=x[2]-x[1])/2
newcost3 = newcost3*100
print('-------')
print('Lowest mismatch Abrasion:', newcost3)
print('Best phi Abrasion:', phi_abr)

# Calculate error
PDF['Art'] = 0.318*PDF['A'] + 0.112*PDF['C'] + 0.207*PDF['F'] + 0.363*PDF['H'] #real values of abrasion

S1 = np.sqrt(PDF['Ami']*PDF['K']).sum()
M1 = np.abs(PDF['Ami']-PDF['K']).sum()/2
L1 = 1-M1

S2 = np.sqrt(PDF['Ami2']*PDF['K']).sum()
M2 = np.abs(PDF['Ami2']-PDF['K']).sum()/2
L2 = 1-M2

S3 = np.sqrt(PDF['Abr']*PDF['K']).sum()
M3 = np.abs(PDF['Abr']-PDF['K']).sum()/2
L3 = 1-M3

S4 = np.sqrt(PDF['Art']*PDF['K']).sum()
M4 = np.abs(PDF['Art']-PDF['K']).sum()/2
L4 = 1-M4

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

rel_error1 = np.trapz(abs(PDF['Art']-PDF['K']).as_matrix().flatten(),
                     dx=x[2]-x[1])/2
rel_error1 = rel_error1*100

# Create distribution
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
for region in objects:
    DIST[region] = stats.rv_discrete(name=region,
                                values=(x, PDF[region].as_matrix().flatten()))
#
# Calculate KS test on two sampled distributions
SIZES['Art'] = SIZES['A'] + SIZES['F'] + SIZES['C'] + SIZES['H']
D1, p1 = 0, 0
D2, p2 = 0, 0
D3, p3 = 0, 0
D4, p4 = 0, 0

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

print('-------')
print('KS 1', kstest1)
print('KS 2', kstest2)
print('KS 3)', kstest3)
print('KS 4', kstest4)
print('-------')
print('KS 1', D1, p1)
print('KS 2', D2, p2)
print('KS 3', D3, p3)
print('KS 4', D4, p4)
print('-------')

# Manual KS Test
def ksprep(data1, data2):
    data1, data2 = map(np.asarray, (data1, data2))
    n1 = data1.shape[0]
    n2 = data2.shape[0]
    n1 = len(data1)
    n2 = len(data2)
    data1 = np.sort(data1)
    data2 = np.sort(data2)
    data_all = np.concatenate([data1,data2])
    cdf1 = np.searchsorted(data1,data_all,side='right')/(1.0*n1)
    cdf2 = (np.searchsorted(data2,data_all,side='right'))/(1.0*n2)

    d = np.max(np.absolute(cdf1-cdf2))
    # Note: d absolute not signed distance
    en = np.sqrt(n1*n2/float(n1+n2))
    prob = stats.distributions.kstwobign.sf(en * d)

    return d, prob

d, en = ksprep(data1, data2)
try:
    prob = stats.distributions.kstwobign.sf(en * d)
except:
    prob = 1.0
print(d, prob)

G1, g1 = 0, 0
G2, g2 = 0, 0
G3, g3 = 0, 0
G4, g4 = 0, 0

for i in range(1000):
    data3, data4 = DIST['Ami'].rvs(size=SIZES['Art']), DIST['K'].rvs(size=SIZES['K'])
    data5, data6 = DIST['Ami2'].rvs(size=SIZES['Art']), DIST['K'].rvs(size=SIZES['K'])
    data7, data8 = DIST['Abr'].rvs(size=SIZES['Art']), DIST['K'].rvs(size=SIZES['K'])
    data1, data2 = DIST['Art'].rvs(size=SIZES['Art']), DIST['K'].rvs(size=SIZES['K'])

    kstest1 = ksprep(data3, data4)
    kstest2 = ksprep(data5, data6)
    kstest3 = ksprep(data7, data8)
    kstest4 = ksprep(data1, data2)
   
    G1 += kstest1[0]
    g1 += kstest1[1]
    G2 += kstest2[0]
    g2 += kstest2[1]
    G3 += kstest3[0]
    g3 += kstest3[1]
    G4 += kstest4[0]
    g4 += kstest4[1]
G1 /= 1000
g1 /= 1000
G2 /= 1000
g2 /= 1000
G3 /= 1000
g3 /= 1000    
G4 /= 1000
g4 /= 1000 

print('KS 1', G1, g1)
print('KS 2', G2, g2)
print('KS 3', G3, g3)
print('KS 4', G4, g4)

print('-------')
#for i in range(1000):
#    # KS for each pop
#    for i, pdf in enumerate(objects2):
#        ni = SIZES[pdf]
#        data1_ = DIST[pdf].rvs(size=ni)
#    
#    for i, pdf in enumerate(objects2):
#        prob_mult = 1
#        d, en = ksprep(data1_, data2)
#        try:
#            prob = stats.distributions.kstwobign.sf(np.sqrt(ni) * d / best_phi[i])
#        except:
#            prob = 1.0
#        prob_mult *= 1-prob
#        print('Region %03s: d=%.2f prob=%f%%' % (pdf, d, 100*prob))
#    print('prob mult=%.5f%%' % (1-prob_mult))

f, ax = plt.subplots()
plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
ax.plot(x, PDF['K'], label='Observed', color = '0.6')
ax.fill_between(x, 0, PDF['K'].as_matrix().flatten(), color = '0.6')
ax.plot(x, PDF['Ami'], 'b',
           label='1')
ax.plot(x, PDF['Ami2'], 'y',
           label='2')
ax.plot(x, PDF['Abr'], 'g',
           label='3')
ax.plot(x, PDF['Art'], 'r',
           label='4')
handles, labels = ax.get_legend_handles_labels()
leg = ax.legend(handles, labels, frameon=True, fancybox=True, fontsize=18)
leg.get_frame().set_edgecolor('k')
ax.set_xlabel('Age [Ma]')
ax.set_ylabel('Relative probability')
ax.set_title('PDPs (X)')
#textstr = 'Mismatch = %.2f%%\nSimilarity = %.2f\nKS test = %.3f' % (rel_error, S, kstest.statistic)
#props = dict(boxstyle='round', facecolor='white', alpha=1.0)
#ax.text(0.7, 0.5, transform=ax.transAxes,
#        verticalalignment='center', bbox=props)
f.savefig(figpath+'pdf_mix_K_E.pdf', bbox_inches='tight')
f.savefig(figpath+'pdf_mix_K_E.png', bbox_inches='tight')

## Create PDF crossplots
x = [PDF['K'].as_matrix().flatten(), PDF['Ami'].as_matrix().flatten(), PDF['Ami2'].as_matrix().flatten(), PDF['Abr'].as_matrix().flatten(), PDF['Art'].as_matrix().flatten()]
#x = {'PDF[E]': (PDF['E'].as_matrix().flatten()), 'PDF[Art]': (PDF['Art'].as_matrix().flatten()), 'PDF[Ami]': (PDF['Ami'].as_matrix().flatten()), 'PDF[Abr]': (PDF['Abr'].as_matrix().flatten()}

gx = 'PDF[K]'
zx = ['PDF[K]', 'PDF[Ami]', 'PDF[Ami2]', 'PDF[Abr]', 'PDF[Art]']
axmin, axmax = min([min(x[0]), min(x[1])]), max([max(x[0]), max(x[1])])
axmin,axmax = axmin-(axmax-axmin)*0.05, axmax+(axmax-axmin)*0.05

# Linear fit
# Don't forget to add ones to fit y = ax+b instead of y = ax
X = sm.add_constant(x[0])
for i in x:
    reg = sm.OLS(i, X).fit()
    sns.set(style='ticks', font_scale=1.5)
    
    # Pretty
    cm = plt.cm.get_cmap('RdYlBu')
    f, ax = plt.subplots()
    ax.plot(x[0], reg.predict(X), 'k', label=r'$R^2 = %.3f$' % reg.rsquared)
    sc = ax.scatter(x[0], i, s=50, alpha=0.6, c=range(4001), facecolors='none', cmap=cm, vmin=0, vmax=4000 )
    ax.set_xlabel(gx)
    ax.set_ylabel('PDF')
    #sns.despine(ax=ax)
    plt.ticklabel_format(style='sci', scilimits=(-2,2))   
    plt.axis('equal')
    ax.set_xlim([axmin, axmax])
    ax.set_ylim([axmin, axmax])
    plt.legend(loc=4)
    cbar=plt.colorbar(sc)
    cbar.solids.set_edgecolor("face")
    plt.draw()
    f.savefig(figpath+'crossplotpdf_K_e.png', bbox_inches='tight')
    f.savefig(figpath+'crossplotpdf_K_e.pdf', bbox_inches='tight')

# Create CDFs
dx = PDF['K'].index[1] - PDF['K'].index[0] 
x = cumtrapz(PDF['K'].as_matrix().flatten(), dx=dx)
y = cumtrapz(PDF['Art'].as_matrix().flatten(), dx=dx)
xh = cumtrapz(PDF['Ami'].as_matrix().flatten(), dx=dx)
xj = cumtrapz(PDF['Ami2'].as_matrix().flatten(), dx=dx)
jj = cumtrapz(PDF['Abr'].as_matrix().flatten(), dx=dx)

axmin, axmax = min([min(x), min(y)]), max([max(x), max(y)])
axmin,axmax = axmin-(axmax-axmin)*0.05, axmax+(axmax-axmin)*0.05
bb = [x, xh, xj, jj, y]
for i in bb:
    X = sm.add_constant(x)
    reg = sm.OLS(i, X).fit()
    
    sns.set(style='ticks', font_scale=1.5)
    f, ax = plt.subplots()
    ax.plot(x, reg.predict(X), 'k', label=r'$R^2 = %.3f$' % reg.rsquared)
    ax.scatter(x, i, s=50, alpha=0.6, edgecolor=prp, facecolors='none')
    ax.set_xlabel('CDF [K]')
    ax.set_ylabel('CDF')
    sns.despine(ax=ax)
    ax.set_xlim([axmin, axmax])
    ax.set_ylim([axmin, axmax])
    plt.legend(loc=4)
    f.savefig(figpath+'crossplotcdf_K_e'+ str(i)+'.png', bbox_inches='tight')
    f.savefig(figpath+'crossplotcdf_K_e.pdf', bbox_inches='tight')

X = sm.add_constant(x)
reg0 = sm.OLS(y, X).fit()
reg1 = sm.OLS(xh, X).fit()
reg2 = sm.OLS(xj, X).fit()
reg3 = sm.OLS(jj, X).fit()

#sns.set(style='ticks', font_scale=1.5)
f, ax = plt.subplots()
ax.plot(x, x, 'k')

ax.scatter(x, y,  s=50, alpha=0.6, edgecolor='r', facecolors='none', label=r'$R^2 = %.3f$' % reg0.rsquared)
ax.scatter(x, xh, s=50, alpha=0.6, edgecolor='b', facecolors='none', label=r'$R^2 = %.3f$' % reg1.rsquared)
ax.scatter(x, xj, s=50, alpha=0.6, edgecolor='y', facecolors='none',label=r'$R^2 = %.3f$' % reg2.rsquared)
ax.scatter(x, jj, s=50, alpha=0.6, edgecolor='g', facecolors='none', label=r'$R^2 = %.3f$' % reg3.rsquared)

ax.set_xlabel('Artificial CDF [K]')
ax.set_ylabel('Predicted CDFs')
ax.set_title('Q-Q plot')
ax.set_xlim([axmin, axmax])
ax.set_ylim([axmin, axmax])
plt.legend(loc=4)
leg.get_frame().set_edgecolor('k')
f.savefig(figpath+'crossplotcdf_G.png', bbox_inches='tight')
f.savefig(figpath+'crossplotcdf_G.pdf', bbox_inches='tight')

x1 = [PDF['K'].as_matrix().flatten(), PDF['Ami'].as_matrix().flatten(), PDF['Ami2'].as_matrix().flatten(), PDF['Abr'].as_matrix().flatten(), PDF['Art'].as_matrix().flatten()]

X = sm.add_constant(x1[0])
reg1 = sm.OLS(x1[1], X).fit()
reg2 = sm.OLS(x1[2], X).fit()
reg3 = sm.OLS(x1[3], X).fit()
reg4 = sm.OLS(x1[4], X).fit()

sns.set(style='ticks', font_scale=1.5)

# Pretty
f, ax = plt.subplots()
ax.plot(x1[0], x1[0], 'k')
ax.scatter(x1[0], x1[1],  s=50, alpha=0.6, edgecolor='b', facecolors='none', label=r'$R^2 = %.3f$' % reg1.rsquared)
ax.scatter(x1[0], x1[2], s=50, alpha=0.6, edgecolor='y', facecolors='none', label=r'$R^2 = %.3f$' % reg2.rsquared)
ax.scatter(x1[0], x1[3], s=50, alpha=0.6, edgecolor='g', facecolors='none',label=r'$R^2 = %.3f$' % reg3.rsquared)
ax.scatter(x1[0], x1[4], s=50, alpha=0.6, edgecolor='r', facecolors='none', label=r'$R^2 = %.3f$' % reg4.rsquared)
ax.set_title('PDF cross-plots')
ax.set_xlabel('PDF')
ax.set_ylabel('PDF')
plt.ticklabel_format(style='sci', scilimits=(-2,2))   
plt.legend(loc=4)
leg.get_frame().set_edgecolor('k')
