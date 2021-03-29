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
# This code expects the original data to be formatted and named properly
################################################################################

# Edit these variables before running
# save_values = 0 # 0-do not save values for statistical tests 1-do save values
make_plots = 1 # 0-do not make plots 1-make plots
make_plot_each_sub = 1 # 0-do not make plots 1-make plots
filter = 1
window = .2

# subjects = [201]
subjects = [201,202,203,204,205,206,207]
# subjects = [203,205]
# number_of_subjects = 3
# min_sub = 101
DIR = "/home/mschlafly/Desktop/Stroke/" #set directory where data is mounted Ola- "/media/ola/Elements/R01prelim" Milli -"Z:"
DT = 0.05
Fs = 1/DT
freq_pendulum = [.5,1,1.5,2.5]
# s103flag = 0 # 0: 10% loading, 1: tabletop

# Label factors as strings for ezANOVA analysis
ind = [0,1,2,3]
frequencylabels = ['0.5Hz','1Hz','1.5Hz','2.5Hz']
# subjects = ['Subject1','Subject2','Subject3','Subject4','Subject5','Subject6','Subject7','Subject8','Subject9','Subject10','Subject11','Subject12','Subject13','Subject14','Subject15','Subject16','Subject17','Subject18']
markerstyles = ['o','D']

# colors = [[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250]]


arms = ['paretic','nonparetic']
arms_label = ['A0','A1']
support_levels = ['0%','35%']
SL_label = ['SL1','SL2']
conditions = ['paretic-0%','paretic-35%','nonparetic-0%','nonparetic-35%']

colors = ['#998ec3', '#f1a340']
linestyles = ['--','-']

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

# if make_plots:
#     xy_xlist = []
#     xy_ylist = []
#     mag_list = []

# if plotting each subject, iterate through each subject plot
num_sub_plots = 1
if make_plot_each_sub:
    num_sub_plots = len(subjects)
for sub_plot in range(0,num_sub_plots):

    # # Iterate through all the files starting with with and without haptic forces
    # compare_groups1 = np.zeros((2,4,number_of_subjects))
    # compare_groups2 = np.zeros((2,4,number_of_subjects))
    # compare_groups3 = np.zeros((len(conditions),4,number_of_subjects))

    for subject_num in range(0,len(subjects)): #iterate though subjects

        if (make_plot_each_sub==0) or (make_plot_each_sub==1 and sub_plot==subject_num):

            if make_plot_each_sub and make_plots:

                # For frequency spectrum plot
                A_Fmag = np.zeros((len(arms),len(support_levels),len(frequencylabels),w_len))
                num_freq = np.zeros((len(arms),len(support_levels),len(frequencylabels)))

                # set up subject scatter plot
                figure_size = (6,3.55)
                fig_sub, ax_sub = plt.subplots(figsize=figure_size,dpi=300)
                xlabel = 'Resonant Frequency of Ball'
                ylabel = 'Fraction of Total Energy'
                title = 'Energy Content of Movement at Different Frequencies'

                # place grid in back
                ax_sub.grid(True, linestyle='-', which='major', axis='y', color='lightgrey',
                               alpha=0.5)
                ax_sub.set_axisbelow(True)

                # Add titles and labels
                plt.xlabel(xlabel,fontname="sans-serif", fontsize=11)
                plt.ylabel(ylabel,fontname="sans-serif", fontsize=11)
                plt.title(title,fontname="sans-serif", fontsize=11,fontweight='bold')
                for label in (ax_sub.get_yticklabels()):
                    label.set_fontsize(8)

                # x-ticks x-axis
                plt.xticks(ind, frequencylabels, fontname="sans-serif", fontsize=10)
                for tick in ax_sub.get_xticklabels():
                    tick.set_rotation(0)

                p = np.zeros(len(arms)+len(support_levels), dtype=object)

            for arm in range(0,len(arms)):

                for SL in range(0,len(support_levels)):

                    # if make_plot_each_sub and make_plots:
                    #     num_freq = np.zeros(4)
                    #     A_Fmag = np.zeros((4,w_len))
                    #     A_Fx = np.zeros((4,w_len))
                    #     A_Fy = np.zeros((4,w_len))

                    arm_SL_list = []

                    for freq in range(0,len(frequencylabels)):

                        energy_at_resonance = np.array([])

                        pathName = DIR + "S" + str(subjects[subject_num]) + "/"
                        Files = []
                        fileNames = os.listdir(pathName)
                        for fileName in fileNames:
                            if fileName.startswith("OutputData"):
                                if "_Freq"+str(freq) in fileName:
                                    if arms_label[arm] in fileName:
                                        if SL_label[SL] in fileName:
                                                # print(fileName)
                                                Files.append(fileName)

                        for i in Files:

                            print(i)
                            file = open(os.path.join(pathName, i), "rU")
                            data = genfromtxt(file,delimiter=',',dtype=float)
                            data = np.delete(data,0,0) # Deletes the first column of column names
                            df = pd.DataFrame({'z':data[:,3],'fx':data[:,8],'fy':data[:,9],'ball_energy':data[:,16]})


                            gravity = 9.81
                            mass = 1.0
                            percent_height = 0.3
                            radius_options = [0.995, 0.249, 0.111, 0.04]
                            R = radius_options[freq]
                            max_energy = mass * gravity * percent_height * R
                            # print('max_energy',max_energy)
                            df = df[df['ball_energy'] > max_energy]

                            table_z = -0.14
                            z_tolerance = table_z + 0.01
                            df = df[df['z'] > z_tolerance]
                            # y_Fx = df['fx'].tolist()
                            # print(hi.tolist())
                            # aslfihaishf
                            # str_split = i.split("_")
                            # trial_num = int(str_split[2].split("Trial")[1])
                            # Use signal directly 1-x position 2- y position 8- x force 9- y force 12-x ball force 13-y ball force
                            # Get amplitudes for force and position signals

                            # y = data[:,8]+data[:,9] # Combine the x and y direction forces
                            # A_Fmag_i,frq = calculate_amplitude(w,y,Fs)
                            # y_Fx = data[:,8]
                            y_Fx = df['fx'].tolist()
                            A_Fx_i,frq = calculate_amplitude(w,y_Fx,Fs)
                            # y_Fy = data[:,9]
                            y_Fy = df['fy'].tolist()
                            A_Fy_i,frq = calculate_amplitude(w,y_Fy,Fs)

                            # Normalize the spectrum so that the energy=1
                            dw = w[1]-w[0]
                            Ax_norm = normalize_spectrum(A_Fx_i,dw)
                            Ay_norm = normalize_spectrum(A_Fy_i,dw)

                            # Add the x and y spectrums together and normalize
                            Amag_norm = normalize_spectrum(Ax_norm+Ay_norm,dw)

                            # Store the signal
                            A_Fmag[arm,SL,freq,:] += Amag_norm
                            num_freq[arm,SL,freq] += 1


                            # A_Fmag[freq,:] += A_Fmag_i
                            # A_Fx[freq,:] += A_Fx_i
                            # A_Fy[freq,:] += A_Fy_i
                            # num_freq[freq] += 1
                            #
                            # dw = w[1]-w[0]
                            # Ax_norm = normalize_spectrum(A_Fx_i,dw)
                            # Ay_norm = normalize_spectrum(A_Fy_i,dw)
                            # Amag_norm = normalize_spectrum(A_Fmag_i,dw)

                            freq_list = []
                            for w_i in range(0,len(w)):
                                if (w[w_i] < freq_pendulum[freq]+window) and (w[w_i] > freq_pendulum[freq]-window):
                                    # freq_list.append(Ax_norm[w_i])
                                    # freq_list.append(Ay_norm[w_i])
                                    freq_list.append(Amag_norm[w_i])

                            energy_at_resonance = np.append(energy_at_resonance,np.sum(np.square(freq_list))*dw)

                            # print(energy_at_resonance)

                        arm_SL_list.append(energy_at_resonance)

                    if make_plot_each_sub and make_plots:

                        # Plot energy@resonance scatterplot
                        data_mean = []
                        data_std = []
                        for freq in range(0,4):
                            data_mean.append(np.mean(arm_SL_list[freq]))
                            data_std.append(np.std(arm_SL_list[freq])/np.sqrt(len(arm_SL_list[freq])))
                            # data_std.append(np.std(arm_SL_list[freq]))

                        # p[group] = ax.errorbar(ind,data_mean,color=colors[group],marker="o")
                        i = arm *2 + SL
                        # print(i)
                        p[i] = ax_sub.errorbar(ind,data_mean,yerr=data_std,color=colors[arm],ls=linestyles[SL],marker="o",ecolor=colors[arm],capsize=5)


                        # print(p)
            if make_plot_each_sub and make_plots:

                # Energy@resonance scatter plot add legend and save
                L = ax_sub.legend(p, conditions, fontsize=10,loc='upper right')
                # L = fig.legend(p, conditions, fontsize=10,loc='upper left',bbox_to_anchor=(0,0.35),bbox_transform=ax.transAxes)#'upper right')
                plt.setp(L.texts, family='sans-serif')
                # fig.subplots_adjust(right=0.7)
                fig_sub.savefig('Plots/IndividualSubjectPlots/'+'S'+str(subjects[subject_num])+'/'+'S'+str(subjects[subject_num])+'.pdf')


                # Make frequency spectrum plots for every ball frequency
                xlabel = ''
                ylabel = 'Frequency Amplitude'
                linestyles_spec = ['-','-','-','-']
                colors_spec = ['#f5793a','#a95aa1','#85c0f9','#0f2080']
                legend = ['Resonant Frequency','paretic-0%','paretic-35%','nonparetic-0%','nonparetic-35%']
                ymin = 0.1
                ymax = 2
                ymax_pend = 1.5
                for freq in range(len(frequencylabels)):
                    mag_list = []
                    for arm in range(0,len(arms)):
                        for SL in range(0,len(support_levels)):
                            A_Fmag[arm,SL,freq,:] /= num_freq[arm,SL,freq]
                            cutoff = 5  # desired cutoff frequency of the filter, Hz
                            A_Fmag[arm,SL,freq,:] = butter_lowpass_filter(A_Fmag[arm,SL,freq,:], cutoff, Fs)
                            A_Fmag[arm,SL,freq,:] = normalize_spectrum(A_Fmag[arm,SL,freq,:],dw)
                            mag_list.append(A_Fmag[arm,SL,freq,:])
                    freq_plot = freq_pendulum#[freq_pendulum[freq]]
                    title = 'Force Frequency Spectrum for the '+frequencylabels[freq]+' Ball'
                    [fig_spec,ax_spec] = mag_spectrum(w,mag_list,freq_plot,title,xlabel,ylabel,legend,linestyles_spec,colors_spec,ymin,ymax,ymax_pend)
                    fig_spec.savefig('Plots/IndividualSubjectPlots/'+'S'+str(subjects[subject_num])+'/'+'S'+str(subjects[subject_num])+'_'+frequencylabels[freq]+'.pdf')

# plt.show()
