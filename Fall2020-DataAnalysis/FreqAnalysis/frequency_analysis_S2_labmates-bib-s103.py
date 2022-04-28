import numpy as np
import scipy as sp
from numpy import fft, genfromtxt
from scipy import fft, arange, stats
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from plot_utils import *
from plot_freq_spectrum import *
import csv
from perform_transform import calculate_amplitude,normalize_spectrum
import pandas as pd
import os

################################################################################
# This program analyses the frequency spectrum for each trial/participant.
# It saves and shows figures and saves matrics into csv file
# "freq-metrics.csv"
# for statistical analyses in R.
# This code expects the original data to be formatted and named properly
################################################################################

subject_num = 103

# Edit these variables before running
save_values = 0 # 0-do not save values for statistical tests 1-do save values
make_plots = 1 # 0-do not make plots 1-make plots
make_plot_each_sub = 1 # 0-do not make plots 1-make plots
haptic_forces_added = 0 # Include adding the haptic forces as an experimental condition
filter = 1
window = .2

number_of_subjects = 3
DIR = "Z:Fall2020Data-round2/" #set directory where data is mounted Ola- "/media/ola/Elements/R01prelim" Milli -"Z:"
DT = 0.05
Fs = 1/DT
freq_pendulum = [.5,1,1.5,2.5]
s103flag = 0 # 0: 10% loading, 1: tabletop

# Label factors as strings for ezANOVA analysis
ind = [0,1,2,3]
frequencylabels = ['0.5Hz','1Hz','1.5Hz','2.5Hz']
markerstyles = ['o','D']


conditions = ['Nonparetic 10% loading','Paretic Tabletop','Paretic 10% Loading']
group_label_read = ['SL2_A1','SL0_A0','SL2_A0']
group_label = ['SL2_A1','SL0_A0','SL2_A0']
colors = [[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250]]

# Create vector of frequencies of interest
nyquist_freq = int(np.floor(Fs/2))
print('Highest frequency evaluated: ', nyquist_freq)
freq_step = 0.1 # set the resolution for the frequency bins

# Fill a vector of frequencies according to these specifications
w = np.zeros(int(np.floor(nyquist_freq/freq_step)-1))
frq = 0 # freq_step
for i in range(len(w)):
    w[i] = frq
    frq = frq + freq_step
print('w: ', w)
w_len = len(w)


# make a numpy array to store data, row=condition col=frequency
data_mat = np.zeros((3,4), dtype=object)

# set up figure
figure_size = (6,4.5) # sets the size of the figure in inches
fig, ax = plt.subplots(figsize=figure_size,dpi=300)
p = np.zeros(len(conditions), dtype=object)

if make_plots:
    xy_xlist = []
    xy_ylist = []
    mag_list = []

for group in range(0,len(conditions)):
# for group in range(0,1):
    num_freq = np.zeros(4)
    A_Fmag = np.zeros((4,w_len))
    A_Fx = np.zeros((4,w_len))
    A_Fy = np.zeros((4,w_len))

    energy_mat = []

    for freq in range(0,4):

        # pathName = DIR + "S0" + str(subject_num) + "/"
        pathName = DIR + "S" + str(subject_num) + "/"
        Files = []
        fileNames = os.listdir(pathName)
        for fileName in fileNames:
            if fileName.startswith("OutputData"):
                if "_Freq"+str(freq) in fileName:
                    if group_label_read[group] in fileName:
                        # print(fileName)
                        Files.append(fileName)

        for i in Files:

            print(i)
            file = open(os.path.join(pathName, i), "rU")
            data = genfromtxt(file,delimiter=',',dtype=float)
            data = np.delete(data,0,0) # Deletes the first column of column names

            str_split = i.split("_")
            trial_num = int(str_split[2].split("Trial")[1])
            if subject_num==103 and trial_num==64:
                continue
            # Use signal directly 1-x position 2- y position 8- x force 9- y force 12-x ball force 13-y ball force
            if 'F1_B1_add'==group_label[group]: # add the haptic forces
                y = data[:,8]+data[:,12]+data[:,9]+data[:,13]
                A_Fmag_i,frq = calculate_amplitude(w,y,Fs)
                y_Fx = data[:,8]+data[:,12]
                A_Fx_i,frq = calculate_amplitude(w,y_Fx,Fs)
                y_Fy = data[:,9]+data[:,13]
                A_Fy_i,frq = calculate_amplitude(w,y_Fy,Fs)
            else:
                # Get amplitudes for force and position signals
                # Combine the x and y direction forces
                y = data[:,8]+data[:,9]
                A_Fmag_i,frq = calculate_amplitude(w,y,Fs)
                y_Fx = data[:,8]
                A_Fx_i,frq = calculate_amplitude(w,y_Fx,Fs)
                y_Fy = data[:,9]
                A_Fy_i,frq = calculate_amplitude(w,y_Fy,Fs)
            A_Fmag[freq,:] += A_Fmag_i
            A_Fx[freq,:] += A_Fx_i
            A_Fy[freq,:] += A_Fy_i
            num_freq[freq] += 1

            dw = w[1]-w[0]
            Ax_norm = normalize_spectrum(A_Fx_i,dw)
            Ay_norm = normalize_spectrum(A_Fy_i,dw)
            Amag_norm = normalize_spectrum(A_Fmag_i,dw)

            # get energy at resonance
            freq_list = []
            for w_i in range(0,len(w)):
                if (w[w_i] < freq_pendulum[freq]+window) and (w[w_i] > freq_pendulum[freq]-window):
                    freq_list.append(Amag_norm[w_i])

            # print(data_mat[group,freq])
            #
            # if type(type(data_mat[group,freq]))=='int':
            #     print('here')
            #     data_mat[group,freq] = np.array(np.sum(np.square(freq_list))*dw)
            #     # energy_mat.append(np.array(np.sum(np.square(freq_list))*dw))
            #     print(type(data_mat[group,freq]))
            # else:
                # print('here')
                # print(np.append(data_mat[group,freq],np.sum(np.square(freq_list))*dw))
                # data_mat[group,freq] = np.append(data_mat[group,freq],np.sum(np.square(freq_list))*dw)
            # if len(data_mat[group,freq]) = 0 = [1,1]
            # data_mat[group,freq].append(np.sum(np.square(freq_list))*dw)
            data_mat[group,freq] = np.append(data_mat[group,freq],np.sum(np.square(freq_list))*dw)

        data_mat[group,freq] = np.delete(data_mat[group,freq],0)
            # print(type(data_mat[group,freq]))
            # # print(energy_mat[freq])
            # try:
            #     if data_mat[group,freq]==0:
            #         data_mat[group,freq] = np.array(np.sum(np.square(freq_list))*dw)
            #         # energy_mat.append(np.array(np.sum(np.square(freq_list))*dw))
            #         print(type(data_mat[group,freq]))
            # except:
            #     print('here')
            #     # print(np.append(data_mat[group,freq],np.sum(np.square(freq_list))*dw))
            #     # data_mat[group,freq] = np.append(data_mat[group,freq],np.sum(np.square(freq_list))*dw)
            #     data_mat[group,freq] = np.append(data_mat[group,freq],np.sum(np.square(freq_list))*dw)

    if make_plots==1:
        data_mean = []
        data_std = []
        for freq in range(0,4):
            data_mean.append(np.mean(data_mat[group,freq]))
            # data_std.append(np.std(data_mat[group,freq])/np.sqrt(data_mat[group,freq].shape[0]))
            data_std.append(np.std(data_mat[group,freq]))
        # print(data_mat[group,:])
        # data_mean.append(=np.mean(data_mat[group,:],axis=1)
        # data_std=np.std(data_mat[group,:])/np.sqrt(data_mat[group,:].shape[0])
        # print(data_mean)
        # print(data_std)

        # ax.errorbar(ind,data_mean)
        p[group] = ax.errorbar(ind,data_mean,color=colors[group],marker="o")
        # p[group] = ax.errorbar(ind,data_mean,yerr=data_std,color=colors[group],marker="o",ecolor=colors[group],capsize=5)




    if make_plots==1:
        # Take average of each signal and re-normalizes
        for freq in range(0,4):
            # if group == 0:
            #     if freq > 0:
            #         A_Fmag[0,:] += A_Fmag[freq,:]
            #         A_Fx[0,:] += A_Fx[freq,:]
            #         A_Fy[0,:] += A_Fy[freq,:]
            #         num_freq[0] += num_freq[freq]
            # else:
            if num_freq[freq]==0:
                print('No data was added for group', group, 'frequency', freq)
                mag_list.append(A_Fmag[freq,:]*0)
                xy_xlist.append(A_Fx[freq,:]*0)
                xy_ylist.append(A_Fy[freq,:]*0)
            else:
                dw = w[1]-w[0]
                A_Fmag[freq,:] /= num_freq[freq]
                A_Fx[freq,:] /= num_freq[freq]
                A_Fy[freq,:] /= num_freq[freq]

                if filter:
                    cutoff = 5  # desired cutoff frequency of the filter, Hz
                    A_Fmag[freq,:] = butter_lowpass_filter(A_Fmag[freq,:], cutoff, Fs)
                    A_Fx[freq,:] = butter_lowpass_filter(A_Fx[freq,:], cutoff, Fs)
                    A_Fy[freq,:] = butter_lowpass_filter(A_Fy[freq,:], cutoff, Fs)

                A_Fmag[freq,:] = normalize_spectrum(A_Fmag[freq,:],dw)
                A_Fx[freq,:] = normalize_spectrum(A_Fx[freq,:],dw)
                A_Fy[freq,:] = normalize_spectrum(A_Fy[freq,:],dw)
                mag_list.append(A_Fmag[freq,:])
                xy_xlist.append(A_Fx[freq,:])
                xy_ylist.append(A_Fy[freq,:])


title = 'Stroke Participant Energy At Resonance'
xlabel = 'Ball\'s Resonant Frequency (Hz)'
ylabel = 'Fraction of Total Energy'

# place grid in back
ax.grid(True, linestyle='-', which='major', axis='y', color='lightgrey',
               alpha=0.5)
ax.set_axisbelow(True)

# Add titles and labels
plt.xlabel(xlabel,fontname="Arial", fontsize=11)
plt.ylabel(ylabel,fontname="Arial", fontsize=11)
plt.title(title,fontname="Arial", fontsize=11,fontweight='bold')
for label in (ax.get_yticklabels()):
    label.set_fontsize(8)

# x-ticks x-axis
plt.xticks(ind, frequencylabels, fontname="Arial", fontsize=10)
for tick in ax.get_xticklabels():
    tick.set_rotation(0)

#figure legend
ax.set_ylim(top=0.4)
L = ax.legend(p, conditions, fontsize=10,loc='upper left')
# L = fig.legend(p, conditions, fontsize=10,loc='upper left',bbox_to_anchor=(0,0.35),bbox_transform=ax.transAxes)#'upper right')
plt.setp(L.texts, family='Arial')
# fig.subplots_adjust(right=0.7)
fig.savefig('Plots/'+'S103_presubmission.pdf')
plt.show()
