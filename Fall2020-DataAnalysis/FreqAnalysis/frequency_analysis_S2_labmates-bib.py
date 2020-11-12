import numpy as np
import scipy as sp
from numpy import fft, genfromtxt
from scipy import fft, arange, stats
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from make_boxplot import make_boxplot, add_stats
import csv
from perform_transform import calculate_amplitude

################################################################################
# This program analyses the frequency spectrum for each trial/participant.
# It saves and shows figures for the BioRob paper and saves the values of the peaks
# around the pendulum's resonant frequency into csv files
# "freq-metrics.csv", "freq-metrics-controls.csv", "freq-metrics-stroke.csv"
# for statistical analyses in R.
# This code expects the original data to be formatted and named properly
################################################################################

# Edit these variables before running
save_values = 1 # 0-do not save 1-save
make_plots = 1 # 0-do not make plots 1-make plots
make_plot_each_sub = 1 # 0-do not make plots 1-make plots
number_of_subjects = 7
DIR = "E:\\Research\\DowntownData\\Fall2020LabData\\"
DT = 0.05
Fs = 1/DT
freq_pendulum = [.5,1,1.5,2]
side_arm = ["R","R","L","L","R","L"]

# Create vector of frequencies of interest
nyquist_freq = int(np.floor(Fs/2))
# nyquist_freq = 40 # override the nyquist frequency (there is little info of value at higher frequencies)
print('Highest frequency evaluated: ', nyquist_freq)
freq_step = 0.1 # set the resolution for the frequcny bins

# Fill a vector of frequencies according to these specifications
w = np.zeros(int(np.floor(nyquist_freq/freq_step)-1))
frq = 0 # freq_step
for i in range(len(w)):
    w[i] = frq
    frq = frq + freq_step
print('w: ', w)
w_len = len(w)

# Label factors as strings for ezANOVA analysis
hapticforces = ['off','on']
frequencylabels = ['same side','cross side']
subjects = ['Subject1','Subject2','Subject3','Subject4','Subject5','Subject6','Subject7','Subject8','Subject9','Subject10','Subject11','Subject12','Subject13','Subject14','Subject15','Subject16','Subject17','Subject18']
colors = ['#601A4A','#EE442F','#63ACBE','#006400',"#E0E0E0"] # for plots
linestyles = ['--','-']
markerstyles = ['o','D']

# Store data for statistical analysis
if save_values==1:
    file = DIR+"freq-metrics.csv"
    columns = ["Subject","HForces","Trial","FmagPeak","DiffPeak","Ratio","FxPeak","FyPeak","PxPeak","PyPeak","RatioFx","RatioFy","FxHigh"]
    with open(file,'w') as csvfile:
        testwriter = csv.writer(csvfile,delimiter=',')
        testwriter.writerow(columns)

freq = 0
# if plotting each subject, iterate through each subject plot
num_iters = 2
if make_plot_each_sub:
    num_iters = number_of_subjects

for sub_plot in range(2,1+num_iters):
    # The average of each amplitude signal is taken. num_0percent-num_50percent records the number of signals added together
    num_freq = np.zeros(2)
    A_Fmag = np.zeros((2,w_len))
    A_Fx = np.zeros((2,w_len))
    A_Fy = np.zeros((2,w_len))
    for subject_num in range(2,1+number_of_subjects): #iterate though subjects

        if (make_plot_each_sub==0) or (make_plot_each_sub==1 and sub_plot==subject_num):
            # Sets the correct subfile
            subname = "S0" + str(subject_num)
            subfile = "\\NailHammerData_S" + str(subject_num)

            for side in range(0,2):
                for k in range(1,14): #iterate through trials
                    try:
                        # Open the trial files
                        trialfile = DIR+subname+subfile+"_Trial"+str(k)+"_Task"+str(side)+"_SL0.csv"
                        data = genfromtxt(trialfile,delimiter=',',dtype=float)
                        data = np.delete(data,0,0) # Deletes the first column of column names
                        print(trialfile)

                        if ((side_arm[subject_num-2]=="R" and side==0) or (side_arm[subject_num-2]=="L" and side==1)):
                            freq = 0
                        else:
                            freq = 1

                        # Get amplitudes for force and position signals
                        # Combine the x and y direction forces
                        y = np.sqrt(np.square(data[:,9])+np.square(data[:,10]))
                        A_Fmag_i,frq = calculate_amplitude(w,y,Fs,True)
                        # print(A_Fmag_i, frq)
                        # print(A_Fmag_i)
                        # isofhaihs
                        # Use signal directly 1-x position 2- y position 8- x force 9- y force
                        y_Fx = data[:,9]
                        A_Fx_i,frq = calculate_amplitude(w,y_Fx,Fs,True)
                        y_Fy = data[:,10]
                        A_Fy_i,frq = calculate_amplitude(w,y_Fy,Fs,True)

                        if make_plots==1:
                            # For combining participants for plots
                            A_Fmag[freq,:] += A_Fmag_i
                            A_Fx[freq,:] += A_Fx_i
                            A_Fy[freq,:] += A_Fy_i
                            num_freq[freq] += 1

                    except:
                        pass

    # if ((make_plot_each_sub==0) or (make_plot_each_sub==1 and plot_num==subject_num)) and make_plots==1:
    if make_plots==1:
        # Take average of each signal
        for freq in range(0,2):
            if num_freq[freq]!=0:
                A_Fmag[freq,:] /= num_freq[freq]
                A_Fx[freq,:] /= num_freq[freq]
                A_Fy[freq,:] /= num_freq[freq]


        # frequency spectrum plot
        wi_0 = 0 # starting index (usually either 0 or 1)
        figure_size = (8,3.55) # inches
        plt.figure(20-sub_plot, figsize=figure_size, dpi=150) #,figsize=(6.0, 4.0), dpi=80)
        # create subplot
        fig_force, ax_force =plt.subplots(nrows=1, ncols=2, sharey='row', squeeze=True, figsize=figure_size, dpi=150)
        fig_force.subplots_adjust(hspace=0.05)
        fig_force.subplots_adjust(wspace=0.05)
        ymin = 0.1
        ymax = 2
        ymax_pend = 1.5
        x1, x2, y1, y2 = 1.2, 2.75, 0.19, 1.75 # specify the limits for inset axis
        zoom = 2.5

        plot_num = 0
        l0 = ax_force[plot_num].plot([0.5, 0.5],[ymin, ymax_pend],linestyle=':',color='k')
        l1 = ax_force[plot_num].plot(w[wi_0:w_len],A_Fx[0,wi_0:w_len],linestyle=linestyles[0],color=colors[0],label=(frequencylabels[0]))
        l2 = ax_force[plot_num].plot(w[wi_0:w_len],A_Fx[1,wi_0:w_len],linestyle=linestyles[0],color=colors[1],label=(frequencylabels[1]))

        # set plot parameters
        ax_force[plot_num].set_xscale('log')
        ax_force[plot_num].set_yscale('log')
        #ax_force[plot_num].set_xlim((10^-1,w[w_len-1]))
        ax_force[plot_num].set_ylim(ymin,ymax)
        ax_force[plot_num].grid(True, color=colors[4])
        ax_force[plot_num].text(.15,.93,'X-Component',horizontalalignment='left',transform=ax_force[plot_num].transAxes, fontsize=10)
        ax_force[plot_num].text(.05,.93,'A.',horizontalalignment='left',transform=ax_force[plot_num].transAxes, fontsize=10, fontweight='bold')


        plot_num = 1
        ax_force[plot_num].plot([0.5, 0.5],[ymin, ymax_pend],linestyle=':',color='k')
        for freq in range(0,2):
            ax_force[plot_num].plot(w[wi_0:w_len],A_Fy[freq,wi_0:w_len],linestyle=linestyles[0],color=colors[freq],label=(frequencylabels[freq]))

        # set plot parameters
        ax_force[plot_num].set_xscale('log')
        ax_force[plot_num].set_yscale('log')
        #ax_force[plot_num].set_xlim((10^-1,w[w_len-1]))
        ax_force[plot_num].set_ylim(ymin,ymax)
        ax_force[plot_num].grid(True, color=colors[4])
        ax_force[plot_num].text(.15,.93,'Y-Component',horizontalalignment='left',transform=ax_force[plot_num].transAxes, fontsize=10)
        ax_force[plot_num].text(.05,.93,'B.',horizontalalignment='left',transform=ax_force[plot_num].transAxes, fontsize=10, fontweight='bold')

        # create titles
        if make_plot_each_sub:
            fig_force.text(0.5, 0.91, 'Force Frequency Spectrum for Subject '+str(sub_plot), ha='center', fontsize=10, fontweight='bold')
        else:
            fig_force.text(0.5, 0.91, 'Force Frequency Spectrums', ha='center', fontsize=10, fontweight='bold')
        # fig_force.suptitle('Force Frequency Spectrums', fontsize=10, fontweight='bold')
        fig_force.text(0.5, 0.01, 'Frequency (Hz)', ha='center', fontsize=10)
        fig_force.text(0.065, 0.5, 'Normalized Frequency Amplitude', va='center', rotation='vertical', fontsize=10)
        # fig_force.legend(loc="upper right", fontsize=9)
        line_labels = ["0.5Hz",
            frequencylabels[0],
            frequencylabels[1]]
        fig_force.legend([l0,l1, l2],labels=line_labels,loc="center right", fontsize=9)

        fig_force.subplots_adjust(right=0.78)

        if make_plot_each_sub:
            fig_force.savefig('freq_forces_S'+str(sub_plot)+'.png')
        else:
            fig_force.savefig('freq_forces.png')

        # # frequency spectrum plot for each frequency
        # freq = 0; group =0
        # figure_size = (8,3.55) # inches
        # plt.figure(21+freq, figsize=figure_size, dpi=150) #,figsize=(6.0, 4.0), dpi=80)
        # if (group==0):
        #     # create subplot
        #     fig_force_freq0, ax_force_freq0 =plt.subplots(nrows=1, ncols=2, sharey='row', squeeze=True, figsize=figure_size, dpi=150)
        #     fig_force_freq0.subplots_adjust(hspace=0.05)
        #     fig_force_freq0.subplots_adjust(wspace=0.05)
        #     ymin = 0.1
        #     ymax = 5
        #     ymax_pend = 1.5
        #     x1, x2, y1, y2 = 1.2, 2.75, 0.19, 1.75 # specify the limits for inset axis
        #     zoom = 2.5
        #
        #     plot_num = 0
        #     l0 = ax_force_freq0[plot_num].plot([freq_pendulum[freq],freq_pendulum[freq]],[ymin, ymax_pend],linestyle=':',color='k')
        #     l1 = ax_force_freq0[plot_num].plot(w[wi_0:w_len],A_Fx[freq,wi_0:w_len],linestyle=linestyles[group],color=colors[freq],label=('forces '+hapticforces[group]))
        #
        #     # set plot parameters
        #     ax_force_freq0[plot_num].set_xscale('log')
        #     ax_force_freq0[plot_num].set_yscale('log')
        #     # ax_force_freq0[plot_num].set_xlim((10^-1,w[w_len-1]))
        #     # ax_force_freq0[plot_num].set_ylim(ymin,ymax)
        #     ax_force_freq0[plot_num].grid(True, color=colors[4])
        #     ax_force_freq0[plot_num].text(.15,.93,'X-Component',horizontalalignment='left',transform=ax_force_freq0[plot_num].transAxes, fontsize=10)
        #     ax_force_freq0[plot_num].text(.05,.93,'A.',horizontalalignment='left',transform=ax_force_freq0[plot_num].transAxes, fontsize=10, fontweight='bold')
        #
        #     plot_num = 1
        #     ax_force_freq0[plot_num].plot([freq_pendulum[freq],freq_pendulum[freq]],[ymin, ymax_pend],linestyle=':',color='k')
        #     ax_force_freq0[plot_num].plot(w[wi_0:w_len],A_Fy[freq,wi_0:w_len],linestyle=linestyles[group],color=colors[freq],label=('forces '+hapticforces[group]))
        #
        #     # set plot parameters
        #     ax_force_freq0[plot_num].set_xscale('log')
        #     ax_force_freq0[plot_num].set_yscale('log')
        #     # ax_force_freq0[plot_num].set_xlim((10^-1,w[w_len-1]))
        #     # ax_force_freq0[plot_num].set_ylim(ymin,ymax)
        #     ax_force_freq0[plot_num].grid(True, color=colors[4])
        #     ax_force_freq0[plot_num].text(.15,.93,'Y-Component',horizontalalignment='left',transform=ax_force_freq0[plot_num].transAxes, fontsize=10)
        #     ax_force_freq0[plot_num].text(.05,.93,'B.',horizontalalignment='left',transform=ax_force_freq0[plot_num].transAxes, fontsize=10, fontweight='bold')
        #
        #
        # freq = 1
        # figure_size = (8,3.55) # inches
        # plt.figure(21+freq, figsize=figure_size, dpi=150) #,figsize=(6.0, 4.0), dpi=80)
        # if (group==0):
        #     # create subplot
        #     fig_force_freq1, ax_force_freq1 =plt.subplots(nrows=1, ncols=2, sharey='row', squeeze=True, figsize=figure_size, dpi=150)
        #     fig_force_freq1.subplots_adjust(hspace=0.05)
        #     fig_force_freq1.subplots_adjust(wspace=0.05)
        #     ymin = 0.1
        #     ymax = 5
        #     ymax_pend = 1.5
        #     x1, x2, y1, y2 = 1.2, 2.75, 0.19, 1.75 # specify the limits for inset axis
        #     zoom = 2.5
        #
        #     plot_num = 0
        #     l0 = ax_force_freq1[plot_num].plot([freq_pendulum[freq],freq_pendulum[freq]],[ymin, ymax_pend],linestyle=':',color='k')
        #     l1 = ax_force_freq1[plot_num].plot(w[wi_0:w_len],A_Fx[freq,wi_0:w_len],linestyle=linestyles[group],color=colors[freq],label=('forces '+hapticforces[group]))
        #
        #     # set plot parameters
        #     ax_force_freq1[plot_num].set_xscale('log')
        #     ax_force_freq1[plot_num].set_yscale('log')
        #     # ax_force_freq1[plot_num].set_xlim((10^-1,w[w_len-1]))
        #     # ax_force_freq1[plot_num].set_ylim(ymin,ymax)
        #     ax_force_freq1[plot_num].grid(True, color=colors[4])
        #     ax_force_freq1[plot_num].text(.15,.93,'X-Component',horizontalalignment='left',transform=ax_force_freq1[plot_num].transAxes, fontsize=10)
        #     ax_force_freq1[plot_num].text(.05,.93,'A.',horizontalalignment='left',transform=ax_force_freq1[plot_num].transAxes, fontsize=10, fontweight='bold')
        #
        #     plot_num = 1
        #     ax_force_freq1[plot_num].plot([freq_pendulum[freq],freq_pendulum[freq]],[ymin, ymax_pend],linestyle=':',color='k')
        #     ax_force_freq1[plot_num].plot(w[wi_0:w_len],A_Fy[freq,wi_0:w_len],linestyle=linestyles[group],color=colors[freq],label=('forces '+hapticforces[group]))
        #
        #     # set plot parameters
        #     ax_force_freq1[plot_num].set_xscale('log')
        #     ax_force_freq1[plot_num].set_yscale('log')
        #     # ax_force_freq1[plot_num].set_xlim((10^-1,w[w_len-1]))
        #     # ax_force_freq1[plot_num].set_ylim(ymin,ymax)
        #     ax_force_freq1[plot_num].grid(True, color=colors[4])
        #     ax_force_freq1[plot_num].text(.15,.93,'Y-Component',horizontalalignment='left',transform=ax_force_freq1[plot_num].transAxes, fontsize=10)
        #     ax_force_freq1[plot_num].text(.05,.93,'B.',horizontalalignment='left',transform=ax_force_freq1[plot_num].transAxes, fontsize=10, fontweight='bold')

plt.show()