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

# subject_num = 7
all_subjects = [1,2,3,4,5,6,7]

# Trial parameters
arm = 0 # 0: paretic; 1: non-paretic
trial_num = 1
support_num = 0 # 0: haptic table (always), 1: 0%, 2: 35%, 3: 50%
freq_num = 0 # 0: 0.5Hz 1: 1Hz 2: 1.5Hz 3: 2.5Hz

# fig, axs = plt.subplots(ncols=3, nrows=2, figsize=(12, 5.5),
#                         constrained_layout=True)
# print(np.array(axs))
# fig.tight_layout()
# np.array(axs).flatten()
# print(axs)
# coordinates = [(0,0), (0,1), (1,0), (1,1), (0,2), (1,2)]

# analysis parameters
DT = 0.05
Fs = 1/DT
freq_pendulum = [.5,1,1.5,2.5]
filter = 1
window = .2

# Create vector of frequencies of interest
nyquist_freq = int(np.floor(Fs/2))
# print('Highest frequency evaluated: ', nyquist_freq)
freq_step = 0.05 # set the resolution for the frequency bins

normalization_cutoff = 8
# Fill a vector of frequencies according to these specifications
w = np.zeros(int(np.floor(nyquist_freq/freq_step)-1))
frq = 0 # freq_step
for i in range(len(w)):
    w[i] = frq
    if frq<normalization_cutoff:
        normalization_cutoff_i = i
    frq = frq + freq_step
w = w[:normalization_cutoff_i]
print('w: ', w)
w_len = len(w)

mag_list = []
A_Fmag = np.zeros(w_len)
for subject_num in all_subjects:
    # Get the location where files are being stored
    path_parent = "E:\\Research\\DowntownWork\\DowntownData\\Fall2020RoundTwo\\" #os.path.dirname(os.path.dirname(os.getcwd()))
    # path_parent = "E:\\Research\\DowntownWork\\DowntownData\\Fall2020RoundTwo\\shaking-ola\\" #os.path.dirname(os.path.dirname(os.getcwd()))
    os.chdir(path_parent)
    pathName = path_parent + "S0" + str(subject_num) + "\\" # + "\\SubjectData\\" #+ "S" + str(subject_num) + "\\"
    fileNames = os.listdir(pathName)

    # Get list of file names that satisfy the desired criterion
    Files = []
    # print(fileNames)
    for fileName in fileNames:
        if fileName.startswith("OutputData"):
            if "_F0_B0" in fileName:
                #if "_A"+str(arm) in fileName:
             #"_SL"+str(support_num) in fileName:
                #print(fileName)
                if all_trials:
                    #Files.append(fileName)
                    # for i in [13,14,15]:
                    # for i in range(13):
                    #if "_Trial"+str(i)+"_" in fileName:
                    Files.append(fileName)
                else:
                    if "_Trial"+str(trial_num) in fileName:
                        Files.append(fileName)
    print(Files)
    if len(Files)<1:
        print('No trials found that match trial parameters')


    # Loop through files
    num_freq = 0

    for i in Files:

        print(i)
        file = open(os.path.join(pathName, i), "rU")
        data = genfromtxt(file,delimiter=',',dtype=float)
        data = np.delete(data,0,0) # Deletes the first column of column names
        df = pd.DataFrame({'z':data[:,3],'fx':data[:,8],'fy':data[:,9],'ball_energy':data[:,11]})
        # df = remove_lowe(df,freq_num)
        # df = remove_notlifted(df)
        y_Fx = df['fx'].tolist()
        A_Fx_i,frq = calculate_amplitude(w,y_Fx,Fs)
        y_Fy = df['fy'].tolist()
        A_Fy_i,frq = calculate_amplitude(w,y_Fy,Fs)

        # plt.plot(frq,A_Fy_i)
        # plt.plot(frq,A_Fx_i)
        # plt.show()

        # Normalize the spectrum so that the energy=1
        dw = w[1]-w[0]
        # Ax_norm = normalize_spectrum(A_Fx_i,dw)
        # Ay_norm = normalize_spectrum(A_Fy_i,dw)

        # Add the x and y spectrums together and normalize
        # Amag_norm = normalize_spectrum(np.sqrt(pow(Ax_norm,2)+pow(Ay_norm,2)),dw)
        # Amag_norm = normalize_spectrum(Ax_norm+Ay_norm,dw)
        Amag_norm = normalize_spectrum(A_Fx_i+A_Fy_i,dw)
        # A_Fmag = A_Fx_i#+A_Fy_i

        # Store the signal
        A_Fmag += Amag_norm
        num_freq += 1

# Make frequency spectrum plot
A_Fmag /= (num_freq*len(all_subjects))
cutoff = 5  # desired cutoff frequency of the filter, Hz
A_Fmag = butter_lowpass_filter(A_Fmag, cutoff, Fs)
A_Fmag = normalize_spectrum(A_Fmag,dw)

# mag_list = [A_Fmag]
mag_list.append(A_Fmag)


###################################################
for subject_num in all_subjects:
    # Get the location where files are being stored
    path_parent = "E:\\Research\\DowntownWork\\DowntownData\\Fall2020RoundTwo\\" #os.path.dirname(os.path.dirname(os.getcwd()))
    # path_parent = "E:\\Research\\DowntownWork\\DowntownData\\Fall2020RoundTwo\\shaking-ola\\" #os.path.dirname(os.path.dirname(os.getcwd()))
    os.chdir(path_parent)
    pathName = path_parent + "S0" + str(subject_num) + "\\" # + "\\SubjectData\\" #+ "S" + str(subject_num) + "\\"
    fileNames = os.listdir(pathName)

    # Get list of file names that satisfy the desired criterion
    Files = []
    # print(fileNames)
    for fileName in fileNames:
        if fileName.startswith("NailHammerData"):
            #if "_Freq"+str(freq_num) in fileName:
                #if "_A"+str(arm) in fileName:
             #"_SL"+str(support_num) in fileName:
                #print(fileName)
            if all_trials:
                #Files.append(fileName)
                for i in [13,14,15]:
                # for i in range(13):
                    if "_Trial"+str(i)+"_" in fileName:
                        Files.append(fileName)
            else:
                if "_Trial"+str(trial_num) in fileName:
                    Files.append(fileName)
    print(Files)
    if len(Files)<1:
        print('No trials found that match trial parameters')
    # analysis parameters
    DT = 0.02
    Fs = 1/DT
    freq_pendulum = [.5,1,1.5,2.5]
    filter = 1
    window = .2

    # Create vector of frequencies of interest
    nyquist_freq = int(np.floor(Fs/2))
    # print('Highest frequency evaluated: ', nyquist_freq)
    freq_step = 0.05 # set the resolution for the frequency bins

    # Fill a vector of frequencies according to these specifications
    # w = np.zeros(int(np.floor(nyquist_freq/freq_step)-1))
    # frq = 0 # freq_step
    # for i in range(len(w)):
    #     w[i] = frq
    #     frq = frq + freq_step
    # # print('w: ', w)
    # w_len = len(w)

    normalization_cutoff = 8
    # Fill a vector of frequencies according to these specifications
    w = np.zeros(int(np.floor(nyquist_freq/freq_step)-1))
    frq = 0 # freq_step
    for i in range(len(w)):
        w[i] = frq
        if frq<normalization_cutoff:
            normalization_cutoff_i = i
        frq = frq + freq_step
    w = w[:normalization_cutoff_i]
    print('w: ', w)
    w_len = len(w)

    # Loop through files
    num_freq = 0
    A_Fmag = np.zeros(w_len)
    for i in Files:

        print(i)
        file = open(os.path.join(pathName, i), "rU")
        data = genfromtxt(file,delimiter=',',dtype=float)
        data = np.delete(data,0,0) # Deletes the first column of column names
        df = pd.DataFrame({'z':data[:,3],'fx':data[:,9],'fy':data[:,10],'ball_energy':data[:,11]})
        # df = remove_lowe(df,freq_num)
        # df = remove_notlifted(df)
        y_Fx = df['fx'].tolist()
        A_Fx_i,frq = calculate_amplitude(w,y_Fx,Fs)
        y_Fy = df['fy'].tolist()
        A_Fy_i,frq = calculate_amplitude(w,y_Fy,Fs)

        # plt.plot(frq,A_Fy_i)
        # plt.plot(frq,A_Fx_i)
        # plt.show()

        # Normalize the spectrum so that the energy=1
        dw = w[1]-w[0]
        # Ax_norm = normalize_spectrum(A_Fx_i,dw)
        # Ay_norm = normalize_spectrum(A_Fy_i,dw)

        # Add the x and y spectrums together and normalize
        # Amag_norm = normalize_spectrum(np.sqrt(pow(Ax_norm,2)+pow(Ay_norm,2)),dw)
        # Amag_norm = normalize_spectrum(Ax_norm+Ay_norm,dw)
        Amag_norm = normalize_spectrum(A_Fx_i+A_Fy_i,dw)
        # A_Fmag = A_Fx_i#+A_Fy_i

        # Store the signal
        A_Fmag += Amag_norm
        num_freq += 1

    # Make frequency spectrum plot
    A_Fmag /= num_freq
    cutoff = 5  # desired cutoff frequency of the filter, Hz
    A_Fmag = butter_lowpass_filter(A_Fmag, cutoff, Fs)
    A_Fmag = normalize_spectrum(A_Fmag,dw)
    # mag_list = [A_Fmag]
    mag_list.append(A_Fmag)


freq_plot = [freq_pendulum[freq_num]]
frequencylabels = ['0.5Hz','1Hz','1.5Hz','2.5Hz']
arms = ['paretic','nonparetic']
support_levels = ['0%','35%']
if all_trials:
    title = '' #'Subject: ' + str(subject_num) #+ ' Arm: ' + arms[arm] + ' SL: ' + support_levels[support_num-1] + ' Freq: ' + frequencylabels[freq_num] + ' Trials: All'
else:
    title = 'Subject: ' + str(subject_num) #+ ' Arm: ' + arms[arm] + ' SL: ' + support_levels[support_num-1] + ' Freq: ' + frequencylabels[freq_num] + ' Trials: ' + str(trial_num)
xlabel = 'Frequency (Hz)'
ylabel = 'Normalized Amplitude'
linestyles_spec = ['-','-','-','-','-','-','-','-']
# colors_spec = ['black','black','black','black','black','black','black']
colors_spec = ['#601A4A','#b2b2b4','#99999b','#7f7f83','#66666a','#4c4c51','#323238','#19191f'] # purple, grey scale
#cccccd
legend = ['','Baseline motion','Subject 1 Max', 'Subject 2 Max','Subject 3 Max', 'Subject 4 Max', 'Subject 5 Max','Subject 6 Max','Subject 7 Max']
ymin = 0
ymax = 1
ymax_pend = 0
[fig_spec,ax_spec] = mag_spectrum(w,mag_list,freq_plot,title,xlabel,ylabel,legend,linestyles_spec,colors_spec,ymin,ymax,ymax_pend)
# axs[coordinates[subject_num-1]] = ax_spec
# plt.savefig('shaking_figure_subject'+str(subject_num)+'.png', bbox_inches="tight")
plt.rcParams["figure.figsize"] = (12,4)

plt.savefig('shaking.svg', bbox_inches="tight")
plt.show()
