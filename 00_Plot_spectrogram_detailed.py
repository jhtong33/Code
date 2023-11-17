import pickle
import glob, os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import matplotlib.dates as mdates
import h5py
from datetime import datetime, date, timedelta
import warnings
warnings.filterwarnings('ignore')


file = '../DataBase/20230414-20230501TaoYuan/wav/PAMGuide_Batch_PSD_Abs_192000ptHannWindow_50pcOlap'

for i, path in enumerate( sorted(glob.glob(f'{file}/7081.23042*.mat'))[1272:1559] ):
    print(path)
    #======= read array 
    arrays = {}
    f = h5py.File(path)
    for k, v in f.items():
        key = k
        arrays[key] = np.array(v)
        
    #======= clean array 
    df_array = arrays[key]
    ra, ca = df_array.shape
    
    time = df_array[0, 1:ca]
    freq = df_array[1:ra, 0]
    df_array = df_array[1:ra, 1:ca]
    
    #============================
    freq_downsample = []
    for f in range(0, ra-1 , 2):
        freq_downsample.append(freq[f])
        temp_row = df_array[f, :]  
        
        templist = []
        for s in range(0, ca-1, 60):
            temp_col = temp_row[s : s + 60]
            pct_95 = np.percentile(temp_col, 95, interpolation='midpoint')
            templist.append(pct_95)
        if f == 0 :
            newarray = np.array(templist)
        else:
            newarray = np.vstack((newarray, np.array(templist)))
    if i == 0 :
        total_array = newarray
    else:
        total_array = np.hstack((total_array, newarray))
    print(total_array.shape)



fig, ax = plt.subplots(1, figsize=(7, 4))
plt.rcParams["font.family"] = "Times New Roman"

vmin = 60
vmax = 120
cmap = cm.jet

freqmin = 1
freqmax = 100

plot_array = total_array[freqmin:freqmax,:]
row, col = plot_array.shape
base = datetime(2023, 4, 24, 10, 0, 0)
alldate = [base + timedelta(minutes = 1*x) for x in range(col)]
print('*'*50)


#========================================================


ax2 = ax.pcolormesh(alldate, np.array(freq_downsample)[freqmin:freqmax], plot_array, cmap=cmap, vmin=vmin, vmax=vmax) 
ax.set_ylabel('Frequency (Hz)', fontsize=15)
ax.set_ylim(freqmin, freqmax)



# cax = fig.add_axes([ax[1].get_position().x1+0.015, ax[1].get_position().y0, 0.01, ax[0].get_position().y0+0.25])
cbar = fig.colorbar(ax2, ax=ax, pad=0.02)
# cbar = plt.colorbar(ax1)
# cbar.ax.set_yticklabels()

cbar.set_label('\nPower Spectral Density\n dB re 1'r'$\mu Pa^2$/Hz', fontsize=12)
cbar.set_ticks(range(vmin,vmax+1,5))


ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))   #to get a tick every 15 minutes
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H'))     #optional formatting 
ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))   #to get a tick every 15 minutes


title = 'Taoyuan_magic_lowfreq_noise_202304'
print(title)
print('-'*70)
plt.title(title)
plt.savefig(f'../Results/{title}.png', dpi = 150)


arr = {'time': alldate, 
       'freq': np.array(freq_downsample),
      '95pct': total_array}

import pickle
pickle.dump( arr, open( f'../DataBase/PAM_pickle/{title}.pkl', 'wb' ) )