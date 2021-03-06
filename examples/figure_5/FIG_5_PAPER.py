# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 20:48:45 2017

@author: Chrystiann
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 20:15:56 2017

@author: s1465002
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

def generate_pdf(filename, sheet, smooth=80, plot=False):
    df = pd.read_excel(filename, sheetname = sheet, skiprows=(0,))

    # Keep only samples whose sample name starts with 504
#    df = df[ df['sample'] != 'SL-13' ]
#    df = df[ df['sample'] != 'xx' ]

    # Create single-column with the true age
    K = df[['Age (Ma)', 'Error (Ma)']]
    K.columns = ['m', 's']
    K.name = 'Sample '+sheet

#    cond = (df['206Pb/238U.1'] + df['206Pb/207Pb.1'])/2 <= 1000
    K['m'] = df['Age (Ma)'] 
    K['s']  = df['Error (Ma)'] 

#    K['m'][~cond] = df[~cond]['206Pb/207Pb.1'] 
#    K['s'][~cond] = df[~cond]['± (Ma)'] 

    # plt.ion()

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
sheets = ['A', 'C', 'F', 'H', 'K']
filename = 'MARSYANDI_SYNTHETIC_1.xlsx'
objects = ('A', 'C', 'K', 'H', 'F')
objects2 = ('A', 'C', 'H', 'F')

PDF = {}
SIZES = {}
for sheet in sheets:
    PDF[sheet], SIZES[sheet] = generate_pdf(filename, sheet, smooth=80, plot=False)
#PDF['C+D'] = (PDF['C'] + PDF['D'])/(PDF['C'] + PDF['D']).sum()
#SIZES['C+D'] = SIZES['C'] + SIZES['D']
#PDF['H+I'] = (PDF['H'] + PDF['I'])/(PDF['H'] + PDF['I']).sum()
#SIZES['H+I'] = SIZES['H'] + SIZES['I']
     
# Experiments with controlling factors against erosion

PDF['Ami'] = 0.525*PDF['A'] + 0.129*PDF['C'] + \
         0.158*PDF['F'] + 0.188*PDF['H'] # abrasion vs erosion

phi_ami2 = [0.654, 0.083, 0.111, 0.152] # Max erosion

PDF['Ami2'] = phi_ami2[0]*PDF['A'] + phi_ami2[1]*PDF['C'] + \
         phi_ami2[2]*PDF['F'] + phi_ami2[3]*PDF['H']

PDF['Abr'] = 0.56*PDF['A'] + 0.1*PDF['C'] + \
         0.14*PDF['F'] + 0.2*PDF['H'] # GSD vs erosion (0.5 for TTS)
         
#PDF['Abr'] = 0.464*PDF['A'] + 0.128*PDF['C'] + 0.173*PDF['F'] + 0.235*PDF['H'] # Max grain-size (0.35 for TTS and 0.15 for the rest)


PDF['Art'] =phi_ami2[0]*PDF['A'] + phi_ami2[1]*PDF['C'] + \
         phi_ami2[2]*PDF['F'] + phi_ami2[3]*PDF['H'] # Max erosion

# Experiments with controlling factors against fertility

PDF['Art3'] = 0.56*PDF['A'] + 0.1*PDF['C'] + \
         0.14*PDF['F'] + 0.2*PDF['H'] # GSD (0.5 for TTS)
         
#PDF['Art3'] = 0.464*PDF['A'] + 0.128*PDF['C'] + 0.173*PDF['F'] + 0.235*PDF['H'] # Max grain-size (0.35 for TTS and 0.15 for the rest)

PDF['Art31'] = 1*PDF['A'] + 0*PDF['C'] + 0*PDF['F'] + 0*PDF['H'] # Max fertility
#PDF['Art31'] = 0.492*PDF['A'] + 0.122*PDF['C'] + 0.164*PDF['F'] + 0.222*PDF['H'] # Max fertility (8.1 for TTS and 3.1 for rest)

phi_abr = [0.654, 0.082, 0.111, 0.153] # Erosion

PDF['Art4'] = phi_abr[0]*PDF['A'] + phi_abr[1]*PDF['C'] + \
         phi_abr[2]*PDF['F'] + phi_abr[3]*PDF['H'] # Erosion
PDF['Art41'] = 1*PDF['A'] + 0*PDF['C'] + 0*PDF['F'] + 0*PDF['H'] # Max fertility
#PDF['Art41'] = 0.492*PDF['A'] + 0.122*PDF['C'] + 0.164*PDF['F'] + 0.222*PDF['H'] # Max fertility (8.1 for TTS and 3.1 for rest)

# Experiments with controlling factors against GSD


#PDF['Art5'] = 0.525*PDF['A'] + 0.129*PDF['C'] + \
#         0.158*PDF['F'] + 0.188*PDF['H'] # abrasion vs gsd
         
PDF['Art5'] = 0.192*PDF['A'] + 0.115*PDF['C'] + \
         0.141*PDF['F'] + 0.552*PDF['H'] # abrasion vs gsd
         
PDF['Art51'] = 0.155*PDF['A'] + 0.1*PDF['C'] + \
         0.135*PDF['F'] + 0.61*PDF['H'] # GSD (0.5 for TTS)
#PDF['Art51'] = 0.464*PDF['A'] + 0.128*PDF['C'] + 0.173*PDF['F'] + 0.235*PDF['H'] # Max grain-size (0.35 for TTS and 0.15 for the rest)

S1 = np.sqrt(PDF['Ami']*PDF['Ami2']).sum()
M1 = np.abs(PDF['Ami']-PDF['Ami2']).sum()/2
L1 = 1-M1

S2 = np.sqrt(PDF['Abr']*PDF['Art']).sum()
M2 = np.abs(PDF['Abr']-PDF['Art']).sum()/2
L2 = 1-M2

S3 = np.sqrt(PDF['Art3']*PDF['Art31']).sum()
M3 = np.abs(PDF['Art3']-PDF['Art31']).sum()/2
L3 = 1-M3

S4 = np.sqrt(PDF['Art4']*PDF['Art41']).sum()
M4 = np.abs(PDF['Art4']-PDF['Art41']).sum()/2
L4 = 1-M4

S5 = np.sqrt(PDF['Art5']*PDF['Art51']).sum()
M5 = np.abs(PDF['Art5']-PDF['Art51']).sum()/2
L5 = 1-M5

print('S1 =', S1)
print('S2 =', S2)
print('S3 =', S3)
print('S4 =', S4)
print('S5 =', S5)


print('-------')
print('M1 =', M1)
print('M2 =', M2)
print('M3 =', M3)
print('M4 =', M4)
print('M5 =', M5)


print('-------')
print('L1 =', L1)
print('L2 =', L2)
print('L3 =', L3)
print('L4 =', L4)
print('L5 =', L5)

x = PDF['Ami'].index

# Create distribution
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
DIST['Art5'] = stats.rv_discrete(name='Art5',
                values=(PDF['Art5'].index, PDF['Art5'].as_matrix().flatten()))
DIST['Art51'] = stats.rv_discrete(name='Art51',
                values=(PDF['Art5'].index, PDF['Art51'].as_matrix().flatten()))
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
D5, p5 = 0, 0

for i in range(1000):
    data1, data2 = DIST['Ami'].rvs(size=SIZES['Art']), DIST['Ami2'].rvs(size=SIZES['Art'])
    data3, data4 = DIST['Abr'].rvs(size=SIZES['Art']), DIST['Art'].rvs(size=SIZES['Art'])
    data5, data6 = DIST['Art3'].rvs(size=SIZES['Art']), DIST['Art31'].rvs(size=SIZES['Art'])
    data7, data8 = DIST['Art4'].rvs(size=SIZES['Art']), DIST['Art41'].rvs(size=SIZES['Art'])
    data9, data10 = DIST['Art5'].rvs(size=SIZES['Art']), DIST['Art51'].rvs(size=SIZES['Art'])
  
    kstest1 = stats.ks_2samp(data1, data2)
    kstest2 = stats.ks_2samp(data3, data4)
    kstest3 = stats.ks_2samp(data5, data6)
    kstest4 = stats.ks_2samp(data7, data8)
    kstest5 = stats.ks_2samp(data9, data10)
   
    D1 += kstest1.statistic
    p1 += kstest1.pvalue
    D2 += kstest2.statistic
    p2 += kstest2.pvalue
    D3 += kstest3.statistic
    p3 += kstest3.pvalue
    D4 += kstest4.statistic
    p4 += kstest4.pvalue
    D5 += kstest5.statistic
    p5 += kstest5.pvalue
    
D1 /= 1000
p1 /= 1000
D2 /= 1000
p2 /= 1000
D3 /= 1000
p3 /= 1000
D4 /= 1000
p4 /= 1000
D5 /= 1000
p5 /= 1000
print('-------')
print('KS 1', kstest1)
print('KS 2', kstest2)
print('KS 3', kstest3)
print('KS 4', kstest4)
print('KS 5', kstest5)

print('-------')
print('KS 1a', D1, p1)
print('KS 2a', D2, p2)
print('KS 3a', D3, p3)
print('KS 4a', D4, p4)
print('KS 5a', D5, p5)

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
G5, g5 = 0, 0
G6, g6 = 0, 0

for i in range(1000):
    data1, data2 = DIST['Ami'].rvs(size=SIZES['Art']), DIST['Ami2'].rvs(size=SIZES['Art'])
    data3, data4 = DIST['Abr'].rvs(size=SIZES['Art']), DIST['Art'].rvs(size=SIZES['Art'])
    data5, data6 = DIST['Art3'].rvs(size=SIZES['Art']), DIST['Art31'].rvs(size=SIZES['Art'])
    data7, data8 = DIST['Art4'].rvs(size=SIZES['Art']), DIST['Art41'].rvs(size=SIZES['Art'])
    data9, data10 = DIST['Art5'].rvs(size=SIZES['Art']), DIST['Art51'].rvs(size=SIZES['Art'])
    data11, data12 = DIST['Art5'].rvs(size=SIZES['Art']), DIST['Art41'].rvs(size=SIZES['Art'])
  
    kstest1 = stats.ks_2samp(data1, data2)
    kstest2 = stats.ks_2samp(data3, data4)
    kstest3 = stats.ks_2samp(data5, data6)
    kstest4 = stats.ks_2samp(data7, data8)
    kstest5 = stats.ks_2samp(data9, data10)
    kstest6 = stats.ks_2samp(data11, data12)
    
    G1 += kstest1[0]
    g1 += kstest1[1]
    G2 += kstest2[0]
    g2 += kstest2[1]
    G3 += kstest3[0]
    g3 += kstest3[1]
    G4 += kstest4[0]
    g4 += kstest4[1]
    G5 += kstest5[0]
    g5 += kstest5[1]
    G6 += kstest6[0]
    g6 += kstest6[1]    
G1 /= 1000
g1 /= 1000
G2 /= 1000
g2 /= 1000
G3 /= 1000
g3 /= 1000
G4 /= 1000
g4 /= 1000
G5 /= 1000
g5 /= 1000
G6 /= 1000
g6 /= 1000
     
print('KS 1', G1, g1)
print('KS 2', G2, g2)
print('KS 3', G3, g3)
print('KS 4', G4, g4)
print('KS 5', G5, g5)
print('KS 6', G6, g6)
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

#f, ax = plt.subplots()
#plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
#ax.plot(x, PDF['Ami'], label='Observed', color = '0.6')
#ax.fill_between(x, 0, PDF['Ami'].as_matrix().flatten(), color = '0.6')
##ax.plot(x, PDF['Ami'], 'b',
##           label='1a')
#ax.plot(x, PDF['Ami2'], 'b',
#           label='1b')
##ax.plot(x, PDF['Abr'], 'g',
##           label='2a')
#ax.plot(x, PDF['Art'], 'g',
#           label='2b')
##ax.plot(x, PDF['Art3'], 'y',
##           label='3a')
#ax.plot(x, PDF['Art31'], 'y',
#           label='3b')
#ax.plot(x, PDF['Art41'], 'c',
#           label='4b')
#handles, labels = ax.get_legend_handles_labels()
#leg = ax.legend(handles, labels, frameon=True, fancybox=True, fontsize=18)
#leg.get_frame().set_edgecolor('k')
#ax.set_xlabel('Age [Ma]')
#ax.set_ylabel('Relative probability')
#ax.set_title('PDPs')
##textstr = 'Mismatch = %.2f%%\nSimilarity = %.2f\nKS test = %.3f' % (rel_error, S, kstest.statistic)
##props = dict(boxstyle='round', facecolor='white', alpha=1.0)
##ax.text(0.7, 0.5, transform=ax.transAxes,
##        verticalalignment='center', bbox=props)
#f.savefig(figpath+'pdf_mix_K_E.pdf', bbox_inches='tight')
#f.savefig(figpath+'pdf_mix_K_E.png', bbox_inches='tight')

## Create PDF crossplots
x = [PDF['Ami'].as_matrix().flatten(), PDF['Ami2'].as_matrix().flatten(), PDF['Abr'].as_matrix().flatten(), PDF['Art'].as_matrix().flatten(), PDF['Art3'].as_matrix().flatten(), PDF['Art31'].as_matrix().flatten(), PDF['Art4'].as_matrix().flatten(), PDF['Art41'].as_matrix().flatten(), PDF['Art5'].as_matrix().flatten(), PDF['Art51'].as_matrix().flatten()]
x1 = [PDF['Ami'].as_matrix().flatten(), PDF['Abr'].as_matrix().flatten(), PDF['Art3'].as_matrix().flatten(), PDF['Art4'].as_matrix().flatten(), PDF['Art5'].as_matrix().flatten()]
x2 = [PDF['Ami2'].as_matrix().flatten(), PDF['Art'].as_matrix().flatten(),PDF['Art31'].as_matrix().flatten(), PDF['Art41'].as_matrix().flatten(), PDF['Art51'].as_matrix().flatten()]


#x = {'PDF[E]': (PDF['E'].as_matrix().flatten()), 'PDF[Art]': (PDF['Art'].as_matrix().flatten()), 'PDF[Ami]': (PDF['Ami'].as_matrix().flatten()), 'PDF[Abr]': (PDF['Abr'].as_matrix().flatten()}

gx1 = ['PDF[Ami]', 'PDF[Abr]', 'PDF[Art3]', 'PDF[Art4]', 'PDF[Art5]']
gx2 = ['PDF[Ami2]','PDF[Art]','PDF[Art31]', 'PDF[Art41]', 'PDF[Art51]']
axmin, axmax = min([min(x[0]), min(x[1])]), max([max(x[0]), max(x[1])])
axmin,axmax = axmin-(axmax-axmin)*0.05, axmax+(axmax-axmin)*0.05

# Linear fit
# Don't forget to add ones to fit y = ax+b instead of y = ax
# X = sm.add_constant(x[2])
#for i in x1:
#    for m in x2:
#        X = sm.add_constant(m)
#        reg = sm.OLS(i, X).fit()
#        sns.set(style='ticks', font_scale=1.5)
#        
#        # Pretty
#        cm = plt.cm.get_cmap('RdYlBu')
#        f, ax = plt.subplots()
#        ax.plot(m, reg.predict(X), 'k', label=r'$R^2 = %.3f$' % reg.rsquared)
#        sc = ax.scatter(m, i, s=50, alpha=0.6, c=range(4001), facecolors='none', cmap=cm, vmin=0, vmax=4000 )
#        ax.set_xlabel(gx1[i])
#        ax.set_ylabel(gx2[i])
#        #sns.despine(ax=ax)
#        plt.ticklabel_format(style='sci', scilimits=(-2,2))   
#        plt.axis('equal')
#        ax.set_xlim([axmin, axmax])
#        ax.set_ylim([axmin, axmax])
#        plt.legend(loc=4)
#        cbar=plt.colorbar(sc)
#        cbar.solids.set_edgecolor("face")
#        plt.draw()
#        f.savefig(figpath+'crossplotpdf_K_e.png', bbox_inches='tight')
#        f.savefig(figpath+'crossplotpdf_K_e.pdf', bbox_inches='tight')

for i,j in zip(x1,x2):
        X = sm.add_constant(j)
        reg = sm.OLS(i, X).fit()
        sns.set(style='ticks', font_scale=1.5)
        
        # Pretty
        cm = plt.cm.get_cmap('RdYlBu')
        f, ax = plt.subplots()
        ax.plot(j, reg.predict(X), 'k', label=r'$R^2 = %.3f$' % reg.rsquared)
        sc = ax.scatter(j, i, s=50, alpha=0.6, c=range(4001), facecolors='none', cmap=cm, vmin=0, vmax=4000 )
        ax.set_xlabel('gx1')
        ax.set_ylabel('gx2')
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
xa = cumtrapz(PDF['Ami'].as_matrix().flatten(), dx=dx)
xb = cumtrapz(PDF['Ami2'].as_matrix().flatten(), dx=dx)
za = cumtrapz(PDF['Abr'].as_matrix().flatten(), dx=dx)
zb = cumtrapz(PDF['Art'].as_matrix().flatten(), dx=dx)
ya = cumtrapz(PDF['Art3'].as_matrix().flatten(), dx=dx)
yb = cumtrapz(PDF['Art31'].as_matrix().flatten(), dx=dx)
ka = cumtrapz(PDF['Art4'].as_matrix().flatten(), dx=dx)
kb = cumtrapz(PDF['Art41'].as_matrix().flatten(), dx=dx)
la = cumtrapz(PDF['Art5'].as_matrix().flatten(), dx=dx)
lb = cumtrapz(PDF['Art51'].as_matrix().flatten(), dx=dx)
axmin, axmax = min([min(xa), min(xb)]), max([max(xa), max(xb)])
axmin,axmax = axmin-(axmax-axmin)*0.05, axmax+(axmax-axmin)*0.05
#bb = [x, y, xh, xj, jj]
#for i in bb:
#    X = sm.add_constant(x)
#    reg = sm.OLS(i, X).fit()
#    
#    sns.set(style='ticks', font_scale=1.5)
#    f, ax = plt.subplots()
#    ax.plot(x, reg.predict(X), 'k', label=r'$R^2 = %.3f$' % reg.rsquared)
#    ax.scatter(x, i, s=50, alpha=0.6, edgecolor=prp, facecolors='none')
#    ax.set_xlabel('CDF [K]')
#    ax.set_ylabel('CDF')
#    sns.despine(ax=ax)
#    ax.set_xlim([axmin, axmax])
#    ax.set_ylim([axmin, axmax])
#    plt.legend(loc=4)
#    f.savefig(figpath+'crossplotcdf_K_e'+ str(i)+'.png', bbox_inches='tight')
#    f.savefig(figpath+'crossplotcdf_K_e.pdf', bbox_inches='tight')

Xa = sm.add_constant(xb)
Xb = sm.add_constant(zb)
Xc = sm.add_constant(yb)
Xd = sm.add_constant(kb)
Xe = sm.add_constant(lb)

reg0 = sm.OLS(xa, Xa).fit()
reg1 = sm.OLS(za, Xb).fit()
reg2 = sm.OLS(ya, Xc).fit()
reg3 = sm.OLS(ka, Xd).fit()
reg4 = sm.OLS(la, Xe).fit()

#sns.set(style='ticks', font_scale=1.5)
f, ax = plt.subplots()
ax.plot(xa, xb, 'k')

ax.scatter(xa, xb,  s=50, alpha=0.6, edgecolor='b', facecolors='none', label=r'$R^2 = %.3f$' % reg0.rsquared)
ax.scatter(za, zb, s=50, alpha=0.6, edgecolor='g', facecolors='none', label=r'$R^2 = %.3f$' % reg1.rsquared)
ax.scatter(ya, yb, s=50, alpha=0.6, edgecolor='y', facecolors='none',label=r'$R^2 = %.3f$' % reg2.rsquared)
ax.scatter(ka, kb, s=50, alpha=0.6, edgecolor='r', facecolors='none',label=r'$R^2 = %.3f$' % reg3.rsquared)
ax.scatter(la, lb, s=50, alpha=0.6, edgecolor='c', facecolors='none',label=r'$R^2 = %.3f$' % reg4.rsquared)

ax.set_xlabel('Abrasion-driven CDFs')
ax.set_ylabel('Erosion-driven CDFs')
ax.set_title('Q-Q plots')
#sns.despine(ax=ax)
ax.set_xlim([axmin, axmax])
ax.set_ylim([axmin, axmax])
plt.legend(loc=4)
#leg.get_frame().set_edgecolor('k')
f.savefig(figpath+'crossplotcdf_G.png', bbox_inches='tight')
f.savefig(figpath+'crossplotcdf_G.pdf', bbox_inches='tight')

# Making PDF cross-plots

x1 = [PDF['Ami'].as_matrix().flatten(), PDF['Abr'].as_matrix().flatten(), PDF['Art3'].as_matrix().flatten(), PDF['Art4'].as_matrix().flatten(), PDF['Art5'].as_matrix().flatten()]
x2 = [PDF['Ami2'].as_matrix().flatten(), PDF['Art'].as_matrix().flatten(),PDF['Art31'].as_matrix().flatten(), PDF['Art41'].as_matrix().flatten(), PDF['Art51'].as_matrix().flatten()]

X1 = sm.add_constant(x1[0])

X2 = sm.add_constant(x2[0])

X3 = sm.add_constant(x2[2])

X4 = sm.add_constant(x2[4])

reg1 = sm.OLS(x1[1], X1).fit() # Ami vs Abr

reg2 = sm.OLS(x1[0], X2).fit() # Ami2 vs Ami
reg3 = sm.OLS(x1[1], X2).fit() # Abr vs Ami

reg4 = sm.OLS(x1[0], X3).fit() # Art31 vs Ami
reg5 = sm.OLS(x1[3], X3).fit() # Art31 vs Art4
reg6 = sm.OLS(x1[1], X3).fit() # Art31 vs Abr

reg7 = sm.OLS(x1[4], X4).fit() # Art5 vs Ami

sns.set(style='ticks', font_scale=1.5)

# Pretty
#cm = plt.cm.get_cmap('RdYlBu')
f, ax = plt.subplots()
#ax.plot(x1[0], x2[0])
#ax.scatter(x1[0], x1[1],  s=50, alpha=0.6, edgecolor='pink', facecolors='none', label=r'$R^2 = %.3f$' % reg1.rsquared)  # Ami vs Abr

ax.scatter(x2[0], x1[0],  s=50, alpha=0.6, edgecolor='b', facecolors='none', label=r'$R^2 = %.3f$' % reg2.rsquared)  # Ami2 vs Ami
ax.scatter(x2[0], x1[1], s=50, alpha=0.6, edgecolor='g', facecolors='none', label=r'$R^2 = %.3f$' % reg3.rsquared) # Abr vs Ami

ax.scatter(x2[2], x1[0], s=50, alpha=0.6, edgecolor='y', facecolors='none',label=r'$R^2 = %.3f$' % reg4.rsquared) # Art31 vs Ami
ax.scatter(x2[2], x1[3], s=50, alpha=0.6, edgecolor='r', facecolors='none',label=r'$R^2 = %.3f$' % reg5.rsquared) # Art31 vs Art4
ax.scatter(x2[2], x1[3], s=50, alpha=0.6, edgecolor='pink', facecolors='none',label=r'$R^2 = %.3f$' % reg6.rsquared) # Art31 vs Art4

ax.scatter(x2[4], x1[0], s=50, alpha=0.6, edgecolor='c', facecolors='none',label=r'$R^2 = %.3f$' % reg7.rsquared) # Art5 vs Ami

#ax.scatter(x[0], i, s=50, alpha=0.6, c=range(4001), facecolors='none', cmap=cm, vmin=0, vmax=4000 )
ax.set_title('PDF cross-plots')
ax.set_xlabel('Abrasion-driven PDFs')
ax.set_ylabel('Erosion-driven PDFs')
#sns.despine(ax=ax)
plt.ticklabel_format(style='sci', scilimits=(-2,2))   
#plt.axis('equal')
#ax.set_xlim([axmin, axmax])
#ax.set_ylim([axmin, axmax])
plt.legend(loc=4)
#leg.get_frame().set_edgecolor('k')

# Three subplots sharing both x/y axes
f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, sharey=False)
#ax1.plot(sample, ero1)
width = 0.10
ax.set_title('PDPs')
ax.set_ylabel('Relative probability')

x = PDF['Ami'].index

ax1.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
ax2.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
ax3.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
ax4.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))

#ax.plot(x, PDF['K'], label='Observed', color = '0.6')
ax1.fill_between(x, 0, PDF['Ami'].as_matrix().flatten(), color = '0.6')
ax2.fill_between(x, 0, PDF['Ami2'].as_matrix().flatten(), color = '0.6')
ax3.fill_between(x, 0, PDF['Art31'].as_matrix().flatten(), color = '0.6')
ax4.fill_between(x, 0, PDF['Art51'].as_matrix().flatten(), color = '0.6')



ax1.plot(x, PDF['Ami'], 'b',
           label='Abrasion')
ax1.plot(x, PDF['Ami'], 'k',
           label='Others')
#ax1.plot(x, PDF['Abr'], 'pink',
#           label='Hillslope GSD')

ax2.plot(x, PDF['Ami'], 'b',
           label='Abrasion')
ax2.plot(x, PDF['Ami2'], 'k',
           label='Erosion')
ax2.plot(x, PDF['Abr'], 'g',
           label='Hillslope GSD')
ax3.plot(x, PDF['Ami'], 'y',
           label='Abrasion')
ax3.plot(x, PDF['Art31'], 'k',
           label='Fertility')
ax3.plot(x, PDF['Art4'], 'r',
           label='Erosion')
ax3.plot(x, PDF['Abr'], 'pink',
           label='Hillslope GSD')
ax4.plot(x, PDF['Art5'], 'c',
           label='Abrasion')
ax4.plot(x, PDF['Art51'], 'k',
           label='Hillslope GSD')


handles, labels = ax.get_legend_handles_labels()
leg = ax.legend(handles, labels, frameon=True, fancybox=True, fontsize=18)
leg.get_frame().set_edgecolor('k')
ax3.set_xlabel('Age [Ma]')
ax3.set_ylabel('Relative probability')
ax1.set_xlim(0, 2000)
# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
#plt.legend()
f.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)


