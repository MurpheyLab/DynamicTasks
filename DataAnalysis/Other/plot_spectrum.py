import numpy as np
import scipy as sp
from numpy import fft, genfromtxt
from scipy import fft, arange, stats
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from plot_utils import *
from plot_freq_utils import *
import csv
from perform_transform import *
import pandas as pd
import os

################################################################################
# This program plots the frequency spectrum for a given trial.
# This code expects the data to be formatted and named properly
################################################################################

all_trials = 1 # 0: only 1 trial 1: all trials for set of conditions

subject_num = 7

# Trial parameters
arm = 0 # 0: paretic; 1: non-paretic
trial_num = 1
support_num = 0 # 0: haptic table (always), 1: 0%, 2: 35%, 3: 50%
freq_num = 0 # 0: 0.5Hz 1: 1Hz 2: 1.5Hz 3: 2.5Hz

# Get the location where files are being stored
path_parent = "E:\\Research\\DowntownWork\\DowntownData\\Fall2020RoundTwo\\shakingdata\\" #os.path.dirname(os.path.dirname(os.getcwd()))
os.chdir(path_parent)
pathName = path_parent + "S0" + str(subject_num) + "\\" # + "\\SubjectData\\" #+ "S" + str(subject_num) + "\\"
fileNames = os.listdir(pathName)
print(fileNames)

# Get list of file names that satisfy the desired criterion
Files = []
for fileName in fileNames:
    if fileName.startswith("OutputData"):
        if "_Freq"+str(freq_num) in fileName:
            #if "_A"+str(arm) in fileName:
            if "_SL"+str(support_num) in fileName:
                print(fileName)
                if all_trials:
                    Files.append(fileName)
                else:
                    if "_Trial"+str(trial_num) in fileName:
                        Files.append(fileName)
print(Files)
if len(Files)<1:
    print('No trials found that match trial parameters')
# analysis parameters
DT = 0.05
Fs = 1/DT
freq_pendulum = [.5,1,1.5,2.5]
filter = 1
window = .2

# Create vector of frequencies of interest
nyquist_freq = int(np.floor(Fs/2))
# print('Highest frequency evaluated: ', nyquist_freq)
freq_step = 0.1 # set the resolution for the frequency bins

# Fill a vector of frequencies according to these specifications
w = np.zeros(int(np.floor(nyquist_freq/freq_step)-1))
frq = 0 # freq_step
for i in range(len(w)):
    w[i] = frq
    frq = frq + freq_step
# print('w: ', w)
w_len = len(w)

# Loop through files
num_freq = 0
A_Fmag = np.zeros(w_len)
for i in Files:

    print(i)
    file = open(os.path.join(pathName, i), "rU")
    data = genfromtxt(file,delimiter=',',dtype=float)
    data = np.delete(data,0,0) # Deletes the first column of column names
    print(data[:,9])
    df = pd.DataFrame({'z':data[:,3],'fx':data[:,9],'fy':data[:,10],'ball_energy':data[:,11]})
    # df = remove_lowe(df,freq_num)
    # df = remove_notlifted(df)
    y_Fx = df['fx'].tolist()
    print("w:", w)
    print("y_Fx:", y_Fx)
    print("Fs:", Fs)
    A_Fx_i,frq = calculate_amplitude(w,y_Fx,Fs)
    y_Fy = df['fy'].tolist()
    A_Fy_i,frq = calculate_amplitude(w,y_Fy,Fs)

    # Normalize the spectrum so that the energy=1
    dw = w[1]-w[0]
    Ax_norm = normalize_spectrum(A_Fx_i,dw)
    Ay_norm = normalize_spectrum(A_Fy_i,dw)

    # Add the x and y spectrums together and normalize
    Amag_norm = normalize_spectrum(Ax_norm+Ay_norm,dw)

    # Store the signal
    A_Fmag += Amag_norm
    num_freq += 1

# Make frequency spectrum plot
A_Fmag /= num_freq
cutoff = 5  # desired cutoff frequency of the filter, Hz
A_Fmag = butter_lowpass_filter(A_Fmag, cutoff, Fs)
A_Fmag = normalize_spectrum(A_Fmag,dw)
mag_list = [A_Fmag]
freq_plot = [freq_pendulum[freq_num]]
frequencylabels = ['0.5Hz','1Hz','1.5Hz','2.5Hz']
arms = ['paretic','nonparetic']
support_levels = ['0%','35%']
if all_trials:
    title = 'Subject: ' + str(subject_num) #+ ' Arm: ' + arms[arm] + ' SL: ' + support_levels[support_num-1] + ' Freq: ' + frequencylabels[freq_num] + ' Trials: All'
else:
    title = 'Subject: ' + str(subject_num) #+ ' Arm: ' + arms[arm] + ' SL: ' + support_levels[support_num-1] + ' Freq: ' + frequencylabels[freq_num] + ' Trials: ' + str(trial_num)
xlabel = ''
ylabel = 'Frequency Amplitude'
linestyles_spec = ['-']
colors_spec = ['black']
legend = []
ymin = 0
ymax = 1
ymax_pend = 0
[fig_spec,ax_spec] = mag_spectrum(w,mag_list,freq_plot,title,xlabel,ylabel,legend,linestyles_spec,colors_spec,ymin,ymax,ymax_pend)
# fig_spec.savefig()
plt.show()
