import csv
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from decimal import Decimal
from scipy.optimize import curve_fit

import scipy.signal as signal

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
#matplotlib.rcParams.update({'font.size': 24})
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)

mup=1

sr=np.linspace(1e-3, 1, num=1e6)
K=1000
mu= mup/sr + mup / sr  *(1-np.exp(-sr*K))

# mup/sr+


major_ticks = np.arange(0, 101, 20)
minor_ticks = np.arange(0, 101, 5)




fig = plt.figure(num=None, figsize=(8, 6), dpi=140, facecolor='w', edgecolor='k')
ax1 = fig.add_subplot(111)
ax1.set_title("Front position", fontsize=26)
ax1.set_ylabel('x($m$)', fontsize=24)
ax1.set_xlabel('time($s$)', fontsize=24)

ax1.grid(color='k', linestyle='-', linewidth=0.2)
ax1.autoscale(enable=True, axis='x', tight=True)
ax1.plot(sr,mu)

leg = ax1.legend()
ax1.set_xlim(0, 1.0)


plt.show()
