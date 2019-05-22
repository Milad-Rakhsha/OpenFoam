import sys
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
# matplotlib.use('Agg')  # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
from decimal import Decimal
from scipy.optimize import curve_fit
import pandas as pd
from collections import OrderedDict

import scipy.signal as signal
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 24})

jop_0 =  ("0.01.csv")

jop_1 =  ("../Jop/0.01_2.csv")
jop_2 =  ("../Jop/0.02_2.csv")
jop_3 =  ("../Jop/0.03_2.csv")
jop_4 =  ("../Jop/0.04_2.csv")
jop_5 =  ("../Jop/0.05_2.csv")




files = OrderedDict([#
#     jop_0: {"method": "Jop", "params":{"name":r'H', "val":0.01} , "shift_t": 0, "lineStyle": '-', "markerEvery": 10},
    (jop_1, {"method": "Jop", "params":{"name":r'H', "val":10} , "shift_t": 0, "lineStyle": '-', "markerEvery": 10}),
    (jop_2, {"method": "Jop", "params":{"name":r'H', "val":20} , "shift_t": 0, "lineStyle": '-', "markerEvery": 10}),
    (jop_3, {"method": "Jop", "params":{"name":r'H', "val":30} , "shift_t": 0, "lineStyle": '-', "markerEvery": 10}),
    (jop_4, {"method": "Jop", "params":{"name":r'H', "val":40} , "shift_t": 0, "lineStyle": '-', "markerEvery": 10}),
    (jop_5, {"method": "Jop", "params":{"name":r'H', "val":50} , "shift_t": 0, "lineStyle": '-', "markerEvery": 10}),
])


major_ticks = np.arange(0, 101, 20)
minor_ticks = np.arange(0, 101, 5)

fig = plt.figure(num=None, figsize=(9, 9),
                 facecolor='w', edgecolor='k')
ax = fig.add_subplot(111)
ax.set_title("Annular Couette ", fontsize=40)
ax.set_ylabel(r'$\omega(r)$', fontsize=40)
ax.set_xlabel(r'$r$ (mm)', fontsize=40)
ax.grid(color='k', linestyle='--', linewidth=0.05)
ax.autoscale(enable=True, axis='x', tight=True)

N=len(files)
print (N)
W=np.zeros(N)
Rc=np.zeros(N)
H=np.zeros(N)

i=0
print (W)
for thisFile in files:
    print (i)
    data = pd.read_csv(thisFile)
    x=data['Points:0']*1000
    y=data['w']
    sr=data['StrainRate']
    eps=5e-3
    x_max=np.max(data[sr>eps]['Points:0'])
    x_min=np.min(data[sr>eps]['Points:0'])
    W[i]=(x_max-x_min)
    Rc_idx=np.argmax(sr)
    Rc[i]=x[Rc_idx]
    H[i]=files[thisFile]["params"]["val"]
    i=i+1
    ax.plot(x, y,
             files[thisFile]["lineStyle"],
             label=r'%s=%2.0f mm' % (
                  files[thisFile]["params"]["name"],files[thisFile]["params"]["val"],
                  ),
             linewidth=1.5, markersize=10,  markevery=files[thisFile]["markerEvery"])

# sub_axes = plt.axes([.2, .6, .2, .2])
# sub_axes.set_title("Shear Band", fontsize=20)
# sub_axes.set_ylabel(r'$R_c$(mm)', fontsize=20)
# sub_axes.set_xlabel(r'$H$ (mm)', fontsize=20)
# sub_axes.plot(H, Rc,'ko', markersize=10)
# sub_axes.set(ylim=[70, 85])
# sub_axes.grid(color='k', linestyle='--', linewidth=0.05)

leg = ax.legend(loc='lower right',
                # bbox_to_anchor=(0.1, -0.55),
                fancybox=True, shadow=True, ncol=1)

# ax.set_facecolor('xkcd:gray')
ax.set_facecolor((0.98, 0.98, 0.98))
# ax.set(xlim=[63, 90])


plt.savefig('Figure_annular_Couette.png',bbox_inches='tight')
plt.show()