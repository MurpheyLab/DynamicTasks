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
DIR = "/media/ola/Data/Research/DowntownData/Fall2020LabData/"
DT = 0.05
Fs = 1/DT
freq_pendulum = [.5,1,1.5,2]


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
frequencylabels = ['0.5Hz','1Hz','1.5Hz','2Hz']
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

# if plotting each subject, iterate through each subject plot
num_iters = 2
if make_plot_each_sub:
    num_iters = number_of_subjects

for sub_plot in range(2,1+num_iters):
    # if make_plot_each_sub:
    #     # For plotting individual subject figures
    #     plt.figure(plot_num)
    #     plt.title('Subject '+str(plot_num)+' Frequency Spectrum')

    # Iterate through all the files starting with with and without haptic forces
    for group in range(len(hapticforces)):

        # The average of each amplitude signal is taken. num_0percent-num_50percent records the number of signals added together
        num_freq = np.zeros(4)
        A_Fmag = np.zeros((4,w_len))
        A_Fx = np.zeros((4,w_len))
        A_Fy = np.zeros((4,w_len))
        for subject_num in range(2,1+number_of_subjects): #iterate though subjects

            if (make_plot_each_sub==0) or (make_plot_each_sub==1 and sub_plot==subject_num):
                # Sets the correct subfile
                subname = "S0" + str(subject_num)
                subfile = "/OutputData_S" + str(subject_num)

                for freq in range(0,4):
                    for k in range(1,50): #iterate through trials
                        try:
                            # if k==12 or k==26 or k==42:
                            #     print('here')
                            #     trialfile = DIR+subname+subfile+"_Task5_SL0_Trial"+str(k)+".csv"
                            #     data = genfromtxt(trialfile,delimiter=',',dtype=float)
                            #     data = np.delete(data,0,0) # Deletes the first column of column names
                            #     print(trialfile)
                            # elif k==7 or k==17 or k==38:
                            #     trialfile = DIR+subname+subfile+"_Task5_SL1_Trial"+str(k)+".csv"
                            #     data = genfromtxt(trialfile,delimiter=',',dtype=float)
                            #     data = np.delete(data,0,0) # Deletes the first column of column names
                            #     print(trialfile)
                            # elif k==5 or k==25 or k==35:
                            #     trialfile = DIR+subname+subfile+"_Task5_SL2_Trial"+str(k)+".csv"
                            #     data = genfromtxt(trialfile,delimiter=',',dtype=float)
                            #     data = np.delete(data,0,0) # Deletes the first column of column names
                            #     print(trialfile)
                            # else:
                            #     hi = hi2

                            # Open the trial files
                            trialfile = DIR+subname+subfile+"_Trial"+str(k)+"_Freq"+str(freq)+"_SL0_F"+str(group)+".csv"
                            # print(trialfile)
                            data = genfromtxt(trialfile,delimiter=',',dtype=float)
                            data = np.delete(data,0,0) # Deletes the first column of column names
                            print(trialfile)

                            # Get amplitudes for force and position signals
                            # Combine the x and y direction forces
                            y = np.sqrt(np.square(data[:,8])+np.square(data[:,9]))
                            A_Fmag_i,frq = calculate_amplitude(w,y,Fs,True)
                            # print(A_Fmag_i)
                            # isofhaihs
                            # Use signal directly 1-x position 2- y position 8- x force 9- y force
                            y_Fx = data[:,8]
                            A_Fx_i,frq = calculate_amplitude(w,y_Fx,Fs,True)
                            y_Fy = data[:,9]
                            A_Fy_i,frq = calculate_amplitude(w,y_Fy,Fs,True)
                            # x_acc = data[delayInd:trial_complete_index,12]
                            # y_Px = y_Fx + x_acc
                            # # y = data[delayInd:trial_complete_index,1] position
                            # A_Px,frq = calculate_amplitude(w,y_Px,Fs,True)
                            # y_acc = data[delayInd:trial_complete_index,12]
                            # y_Py = y_Fy + y_acc
                            # # y = data[delayInd:trial_complete_index,2] position
                            # A_Py,frq = calculate_amplitude(w,y_Py,Fs,True)
                            # y_Fmag_net = np.sqrt(np.square(y_Px)+np.square(y_Py))
                            # A_Fmag_net,frq = calculate_amplitude(w,y_Fmag_net,Fs,True)
                            # A_Dx = A_Fx - A_Px
                            # A_Dy = A_Fy - A_Py
                            # A_Fmag_diff = A_Fmag - A_Fmag_net

                            if make_plots==1:
                                # For combining participants for plots
                                A_Fmag[freq,:] += A_Fmag_i
                                A_Fx[freq,:] += A_Fx_i
                                A_Fy[freq,:] += A_Fy_i
                                num_freq[freq] += 1

                        except:
                            pass
                    # if make_plot_each_sub:
                    #     # For individual subject plot
                    #     A_support /= num_supportlevel
                    #     plt.plot(w[:w_len],A_support[:w_len],linestyle=linestyles[group],color=colors[j],label=supportlevels[j]) # plotting the spectrum (start at )
                    #     plt.legend(loc="upper right")
                    #     plt.xlabel('Freq (Hz)')
                    #     plt.ylabel('Amplitude')
                    #     plt.xscale('log')
                    #     plt.yscale('log')
                    #     plt.savefig('sub'+str(subject_num)+'_freq.png')
        # if ((make_plot_each_sub==0) or (make_plot_each_sub==1 and plot_num==subject_num)) and make_plots==1:
        if make_plots==1:
            # Take average of each signal
            for freq in range(0,4):
                if num_freq[freq]!=0:
                    A_Fmag[freq,:] /= num_freq[freq]
                    A_Fx[freq,:] /= num_freq[freq]
                    A_Fy[freq,:] /= num_freq[freq]


            # Figures with both haptic force levels

            # frequency spectrum plot
            wi_0 = 0 # starting index (usually either 0 or 1)
            figure_size = (8,3.55) # inches
            plt.figure(20-sub_plot, figsize=figure_size, dpi=150) #,figsize=(6.0, 4.0), dpi=80)
            if (group==0):
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
                l1 = ax_force[plot_num].plot(w[wi_0:w_len],A_Fx[0,wi_0:w_len],linestyle=linestyles[group],color=colors[0],label=('forces '+hapticforces[group]+' '+frequencylabels[0]))
                l2 = ax_force[plot_num].plot(w[wi_0:w_len],A_Fx[1,wi_0:w_len],linestyle=linestyles[group],color=colors[1],label=('forces '+hapticforces[group]+' '+frequencylabels[1]))
                l3 = ax_force[plot_num].plot(w[wi_0:w_len],A_Fx[2,wi_0:w_len],linestyle=linestyles[group],color=colors[2],label=('forces '+hapticforces[group]+' '+frequencylabels[2]))
                l4 = ax_force[plot_num].plot(w[wi_0:w_len],A_Fx[3,wi_0:w_len],linestyle=linestyles[group],color=colors[3],label=('forces '+hapticforces[group]+' '+frequencylabels[3]))

                # set plot parameters
                ax_force[plot_num].set_xscale('log')
                ax_force[plot_num].set_yscale('log')
                ax_force[plot_num].set_xlim((10^-1,w[w_len-1]))
                ax_force[plot_num].set_ylim(ymin,ymax)
                ax_force[plot_num].grid(True, color=colors[4])
                ax_force[plot_num].text(.15,.93,'X-Component',horizontalalignment='left',transform=ax_force[plot_num].transAxes, fontsize=10)
                ax_force[plot_num].text(.05,.93,'A.',horizontalalignment='left',transform=ax_force[plot_num].transAxes, fontsize=10, fontweight='bold')
                # # create inset axes
                # axins_for0 = zoomed_inset_axes(ax_force[plot_num], zoom, loc='lower left')#,bbox_to_anchor=(0.5, 8, 4, 4)) # bbox_to_anchor=(10,8)) # zoom-factor: 2.5, location: upper-left
                # axins_for0.set_xlim(x1, x2) # apply the x-limits
                # axins_for0.set_ylim(y1, y2) # apply the y-limits
                # plt.tick_params(
                #     axis='both',
                #     which='both',
                #     bottom=False,
                #     left=False,
                #     right=False,
                #     top=False,
                #     labelleft=False,
                #     labelbottom=False)
                # axins_for0.set_xscale('log')
                # axins_for0.set_yscale('log')
                # axins_for0.grid(True, color=colors[3], which='both')
                # mark_inset(ax_force[plot_num], axins_for0, loc1=2, loc2=4, fc="none", ec="0.5",edgecolor='#000000', facecolor='#000000', color='#000000')
                # axins_for0.plot([freq_pendulum, freq_pendulum],[ymin, ymax],linestyle=':',color='k')

                plot_num = 1
                ax_force[plot_num].plot([0.5, 0.5],[ymin, ymax_pend],linestyle=':',color='k')
                for freq in range(0,4):
                    ax_force[plot_num].plot(w[wi_0:w_len],A_Fy[freq,wi_0:w_len],linestyle=linestyles[group],color=colors[freq],label=('forces '+hapticforces[group]+' '+frequencylabels[freq]))

                # set plot parameters
                ax_force[plot_num].set_xscale('log')
                ax_force[plot_num].set_yscale('log')
                ax_force[plot_num].set_xlim((10^-1,w[w_len-1]))
                ax_force[plot_num].set_ylim(ymin,ymax)
                ax_force[plot_num].grid(True, color=colors[4])
                ax_force[plot_num].text(.15,.93,'Y-Component',horizontalalignment='left',transform=ax_force[plot_num].transAxes, fontsize=10)
                ax_force[plot_num].text(.05,.93,'B.',horizontalalignment='left',transform=ax_force[plot_num].transAxes, fontsize=10, fontweight='bold')
                # # create inset axes
                # axins_for1 = zoomed_inset_axes(ax_force[plot_num], zoom, loc='lower left')#,bbox_to_anchor=(0.5, 8, 4, 4)) # bbox_to_anchor=(10,8)) # zoom-factor: 2.5, location: upper-left
                # axins_for1.set_xlim(x1, x2) # apply the x-limits
                # axins_for1.set_ylim(y1, y2) # apply the y-limits
                # plt.tick_params(
                #     axis='both',
                #     which='both',
                #     bottom=False,
                #     left=False,
                #     right=False,
                #     top=False,
                #     labelleft=False,
                #     labelbottom=False)
                # axins_for1.set_xscale('log')
                # axins_for1.set_yscale('log')
                # axins_for1.grid(True, color=colors[3], which='both')
                # mark_inset(ax_force[plot_num], axins_for1, loc1=2, loc2=4, fc="none", ec="0.5",edgecolor='#000000', facecolor='#000000', color='#000000')
                # axins_for1.plot([freq_pendulum, freq_pendulum],[ymin, ymax],linestyle=':',color='k')

            # axins_for0.plot(w[wi_0:w_len],A_Fx_0percent[wi_0:w_len],linestyle=linestyles[group],color=colors[0],label=(group_name[group]+' '+supportlevels[0]))
            # axins_for0.plot(w[wi_0:w_len],A_Fx_20percent[wi_0:w_len],linestyle=linestyles[group],color=colors[1],label=(group_name[group]+' '+supportlevels[1]))
            # axins_for0.plot(w[wi_0:w_len],A_Fx_50percent[wi_0:w_len],linestyle=linestyles[group],color=colors[2],label=(group_name[group]+' '+supportlevels[2]))
            # axins_for1.plot(w[wi_0:w_len],A_Fy_0percent[wi_0:w_len],linestyle=linestyles[group],color=colors[0],label=(group_name[group]+' '+supportlevels[0]))
            # axins_for1.plot(w[wi_0:w_len],A_Fy_20percent[wi_0:w_len],linestyle=linestyles[group],color=colors[1],label=(group_name[group]+' '+supportlevels[1]))
            # axins_for1.plot(w[wi_0:w_len],A_Fy_50percent[wi_0:w_len],linestyle=linestyles[group],color=colors[2],label=(group_name[group]+' '+supportlevels[2]))

            if (group==1):
                plot_num = 0
                l5 = ax_force[plot_num].plot(w[wi_0:w_len],A_Fx[0,wi_0:w_len],linestyle=linestyles[group],color=colors[0],label=('forces '+hapticforces[group]+' '+frequencylabels[0]))
                l6 = ax_force[plot_num].plot(w[wi_0:w_len],A_Fx[1,wi_0:w_len],linestyle=linestyles[group],color=colors[1],label=('forces '+hapticforces[group]+' '+frequencylabels[1]))
                l7 = ax_force[plot_num].plot(w[wi_0:w_len],A_Fx[2,wi_0:w_len],linestyle=linestyles[group],color=colors[2],label=('forces '+hapticforces[group]+' '+frequencylabels[2]))
                l8 = ax_force[plot_num].plot(w[wi_0:w_len],A_Fx[3,wi_0:w_len],linestyle=linestyles[group],color=colors[3],label=('forces '+hapticforces[group]+' '+frequencylabels[3]))
                ax_force[plot_num].plot([1, 1],[ymin, ymax_pend],linestyle=':',color='k')
                ax_force[plot_num].plot([1.5, 1.5],[ymin, ymax_pend],linestyle=':',color='k')
                ax_force[plot_num].plot([2, 2],[ymin, ymax_pend],linestyle=':',color='k')

                # set plot parameters

                for label in (ax_force[plot_num].get_xticklabels() + ax_force[plot_num].get_yticklabels()):
                    label.set_fontsize(8)

                plot_num = 1
                ax_force[plot_num].plot([1, 1],[ymin, ymax_pend],linestyle=':',color='k')
                ax_force[plot_num].plot([1.5, 1.5],[ymin, ymax_pend],linestyle=':',color='k')
                ax_force[plot_num].plot([2, 2],[ymin, ymax_pend],linestyle=':',color='k')
                for freq in range(0,4):
                    ax_force[plot_num].plot(w[wi_0:w_len],A_Fy[freq,wi_0:w_len],linestyle=linestyles[group],color=colors[freq],label=('forces '+hapticforces[group]+' '+frequencylabels[freq]))
                for label in (ax_force[plot_num].get_xticklabels() + ax_force[plot_num].get_yticklabels()):
                    label.set_fontsize(8)

                # create titles
                if make_plot_each_sub:
                    fig_force.text(0.5, 0.91, 'Force Frequency Spectrum for Subject '+str(sub_plot), ha='center', fontsize=10, fontweight='bold')
                else:
                    fig_force.text(0.5, 0.91, 'Force Frequency Spectrums', ha='center', fontsize=10, fontweight='bold')
                # fig_force.suptitle('Force Frequency Spectrums', fontsize=10, fontweight='bold')
                fig_force.text(0.5, 0.01, 'Frequency (Hz)', ha='center', fontsize=10)
                fig_force.text(0.065, 0.5, 'Normalized Frequency Amplitude', va='center', rotation='vertical', fontsize=10)
                # fig_force.legend(loc="upper right", fontsize=9)
                line_labels = ["Ball's resonant\nfrequency",
                    'forces '+hapticforces[0]+' '+frequencylabels[0],
                    'forces '+hapticforces[0]+' '+frequencylabels[1],
                    'forces '+hapticforces[0]+' '+frequencylabels[2],
                    'forces '+hapticforces[0]+' '+frequencylabels[3],
                    'forces '+hapticforces[1]+' '+frequencylabels[0],
                    'forces '+hapticforces[1]+' '+frequencylabels[1],
                    'forces '+hapticforces[1]+' '+frequencylabels[2],
                    'forces '+hapticforces[1]+' '+frequencylabels[3]]
                #["Ball's resonant\nfrequency",group_name[0]+' '+supportlevels[0],group_name[0]+' '+supportlevels[1],group_name[0]+' '+supportlevels[2],group_name[1]+' '+supportlevels[0],group_name[1]+' '+supportlevels[1],group_name[1]+' '+supportlevels[2]]
                fig_force.legend([l0,l1, l2, l3, l4, l5, l6, l7, l8],labels=line_labels,loc="center right", fontsize=9)
                # fig.legend([l1, l2, l3, l4, l5, l6],labels=line_labels,loc="center right", fontsize=9)
                fig_force.subplots_adjust(right=0.78)
                # label_x = 8
                # label_y = 0.0007
                # arrow_x = freq_pendulum*2
                # arrow_y = 0.003
                #
                # label_x = 18
                # label_y = 0.1
                # arrow_x = freq_pendulum*2
                # arrow_y = 0.6
                # arrow_properties = dict(
                #     facecolor="black", width=0.5,
                #     headwidth=4, shrink=0.1)
                # ax_force.annotate(
                #     "Pendulum's\nresonant\nfrequency", xy=(arrow_x, arrow_y),
                #     xytext=(label_x, label_y),
                #     arrowprops=arrow_properties, ha='center')

                if make_plot_each_sub:
                    fig_force.savefig('freq_forces_S'+str(sub_plot)+'.png')
                else:
                    fig_force.savefig('freq_forces.png')


            # frequency spectrum plot for each frequency
            freq = 0
            figure_size = (8,3.55) # inches
            plt.figure(21+freq, figsize=figure_size, dpi=150) #,figsize=(6.0, 4.0), dpi=80)
            if (group==0):
                # create subplot
                fig_force_freq0, ax_force_freq0 =plt.subplots(nrows=1, ncols=2, sharey='row', squeeze=True, figsize=figure_size, dpi=150)
                fig_force_freq0.subplots_adjust(hspace=0.05)
                fig_force_freq0.subplots_adjust(wspace=0.05)
                ymin = 0.1
                ymax = 5
                ymax_pend = 1.5
                x1, x2, y1, y2 = 1.2, 2.75, 0.19, 1.75 # specify the limits for inset axis
                zoom = 2.5

                plot_num = 0
                l0 = ax_force_freq0[plot_num].plot([freq_pendulum[freq],freq_pendulum[freq]],[ymin, ymax_pend],linestyle=':',color='k')
                l1 = ax_force_freq0[plot_num].plot(w[wi_0:w_len],A_Fx[freq,wi_0:w_len],linestyle=linestyles[group],color=colors[freq],label=('forces '+hapticforces[group]))

                # set plot parameters
                ax_force_freq0[plot_num].set_xscale('log')
                ax_force_freq0[plot_num].set_yscale('log')
                # ax_force_freq0[plot_num].set_xlim((10^-1,w[w_len-1]))
                # ax_force_freq0[plot_num].set_ylim(ymin,ymax)
                ax_force_freq0[plot_num].grid(True, color=colors[4])
                ax_force_freq0[plot_num].text(.15,.93,'X-Component',horizontalalignment='left',transform=ax_force_freq0[plot_num].transAxes, fontsize=10)
                ax_force_freq0[plot_num].text(.05,.93,'A.',horizontalalignment='left',transform=ax_force_freq0[plot_num].transAxes, fontsize=10, fontweight='bold')

                plot_num = 1
                ax_force_freq0[plot_num].plot([freq_pendulum[freq],freq_pendulum[freq]],[ymin, ymax_pend],linestyle=':',color='k')
                ax_force_freq0[plot_num].plot(w[wi_0:w_len],A_Fy[freq,wi_0:w_len],linestyle=linestyles[group],color=colors[freq],label=('forces '+hapticforces[group]))

                # set plot parameters
                ax_force_freq0[plot_num].set_xscale('log')
                ax_force_freq0[plot_num].set_yscale('log')
                # ax_force_freq0[plot_num].set_xlim((10^-1,w[w_len-1]))
                # ax_force_freq0[plot_num].set_ylim(ymin,ymax)
                ax_force_freq0[plot_num].grid(True, color=colors[4])
                ax_force_freq0[plot_num].text(.15,.93,'Y-Component',horizontalalignment='left',transform=ax_force_freq0[plot_num].transAxes, fontsize=10)
                ax_force_freq0[plot_num].text(.05,.93,'B.',horizontalalignment='left',transform=ax_force_freq0[plot_num].transAxes, fontsize=10, fontweight='bold')

            if (group==1):
                plot_num = 0
                l5 = ax_force_freq0[plot_num].plot(w[wi_0:w_len],A_Fx[freq,wi_0:w_len],linestyle=linestyles[group],color=colors[0],label=('forces '+hapticforces[group]))
                for label in (ax_force_freq0[plot_num].get_xticklabels() + ax_force_freq0[plot_num].get_yticklabels()):
                    label.set_fontsize(8)

                plot_num = 1
                ax_force_freq0[plot_num].plot(w[wi_0:w_len],A_Fy[freq,wi_0:w_len],linestyle=linestyles[group],color=colors[freq],label=('forces '+hapticforces[group]+' '+frequencylabels[freq]))
                for label in (ax_force_freq0[plot_num].get_xticklabels() + ax_force_freq0[plot_num].get_yticklabels()):
                    label.set_fontsize(8)

                # create titles
                fig_force_freq0.text(0.5, 0.91, 'Force Frequency Spectrums for 0.5Hz Ball', ha='center', fontsize=10, fontweight='bold')
                # freq.suptitle('Force Frequency Spectrums', fontsize=10, fontweight='bold')
                fig_force_freq0.text(0.5, 0.01, 'Frequency (Hz)', ha='center', fontsize=10)
                fig_force_freq0.text(0.065, 0.5, 'Normalized Frequency Amplitude', va='center', rotation='vertical', fontsize=10)
                # freq.legend(loc="upper right", fontsize=9)
                line_labels = ["Ball's resonant\nfrequency",
                    'forces '+hapticforces[0],
                    'forces '+hapticforces[1]]
                fig_force_freq0.legend([l0, l1, l2],labels=line_labels,loc="center right", fontsize=9)
                fig_force_freq0.subplots_adjust(right=0.78)

                fig_force_freq0.savefig('freq_forces'+str(freq)+'.png')
            freq = 1
            figure_size = (8,3.55) # inches
            plt.figure(21+freq, figsize=figure_size, dpi=150) #,figsize=(6.0, 4.0), dpi=80)
            if (group==0):
                # create subplot
                fig_force_freq1, ax_force_freq1 =plt.subplots(nrows=1, ncols=2, sharey='row', squeeze=True, figsize=figure_size, dpi=150)
                fig_force_freq1.subplots_adjust(hspace=0.05)
                fig_force_freq1.subplots_adjust(wspace=0.05)
                ymin = 0.1
                ymax = 5
                ymax_pend = 1.5
                x1, x2, y1, y2 = 1.2, 2.75, 0.19, 1.75 # specify the limits for inset axis
                zoom = 2.5

                plot_num = 0
                l0 = ax_force_freq1[plot_num].plot([freq_pendulum[freq],freq_pendulum[freq]],[ymin, ymax_pend],linestyle=':',color='k')
                l1 = ax_force_freq1[plot_num].plot(w[wi_0:w_len],A_Fx[freq,wi_0:w_len],linestyle=linestyles[group],color=colors[freq],label=('forces '+hapticforces[group]))

                # set plot parameters
                ax_force_freq1[plot_num].set_xscale('log')
                ax_force_freq1[plot_num].set_yscale('log')
                # ax_force_freq1[plot_num].set_xlim((10^-1,w[w_len-1]))
                # ax_force_freq1[plot_num].set_ylim(ymin,ymax)
                ax_force_freq1[plot_num].grid(True, color=colors[4])
                ax_force_freq1[plot_num].text(.15,.93,'X-Component',horizontalalignment='left',transform=ax_force_freq1[plot_num].transAxes, fontsize=10)
                ax_force_freq1[plot_num].text(.05,.93,'A.',horizontalalignment='left',transform=ax_force_freq1[plot_num].transAxes, fontsize=10, fontweight='bold')

                plot_num = 1
                ax_force_freq1[plot_num].plot([freq_pendulum[freq],freq_pendulum[freq]],[ymin, ymax_pend],linestyle=':',color='k')
                ax_force_freq1[plot_num].plot(w[wi_0:w_len],A_Fy[freq,wi_0:w_len],linestyle=linestyles[group],color=colors[freq],label=('forces '+hapticforces[group]))

                # set plot parameters
                ax_force_freq1[plot_num].set_xscale('log')
                ax_force_freq1[plot_num].set_yscale('log')
                # ax_force_freq1[plot_num].set_xlim((10^-1,w[w_len-1]))
                # ax_force_freq1[plot_num].set_ylim(ymin,ymax)
                ax_force_freq1[plot_num].grid(True, color=colors[4])
                ax_force_freq1[plot_num].text(.15,.93,'Y-Component',horizontalalignment='left',transform=ax_force_freq1[plot_num].transAxes, fontsize=10)
                ax_force_freq1[plot_num].text(.05,.93,'B.',horizontalalignment='left',transform=ax_force_freq1[plot_num].transAxes, fontsize=10, fontweight='bold')

            if (group==1):
                plot_num = 0
                l5 = ax_force_freq1[plot_num].plot(w[wi_0:w_len],A_Fx[freq,wi_0:w_len],linestyle=linestyles[group],color=colors[freq],label=('forces '+hapticforces[group]))
                for label in (ax_force_freq1[plot_num].get_xticklabels() + ax_force_freq1[plot_num].get_yticklabels()):
                    label.set_fontsize(8)

                plot_num = 1
                ax_force_freq1[plot_num].plot(w[wi_0:w_len],A_Fy[freq,wi_0:w_len],linestyle=linestyles[group],color=colors[freq],label=('forces '+hapticforces[group]+' '+frequencylabels[freq]))
                for label in (ax_force_freq1[plot_num].get_xticklabels() + ax_force_freq1[plot_num].get_yticklabels()):
                    label.set_fontsize(8)

                # create titles
                fig_force_freq1.text(0.5, 0.91, 'Force Frequency Spectrums for 1Hz Ball', ha='center', fontsize=10, fontweight='bold')
                # freq.suptitle('Force Frequency Spectrums', fontsize=10, fontweight='bold')
                fig_force_freq1.text(0.5, 0.01, 'Frequency (Hz)', ha='center', fontsize=10)
                fig_force_freq1.text(0.065, 0.5, 'Normalized Frequency Amplitude', va='center', rotation='vertical', fontsize=10)
                # freq.legend(loc="upper right", fontsize=9)
                line_labels = ["Ball's resonant\nfrequency",
                    'forces '+hapticforces[0],
                    'forces '+hapticforces[1]]
                fig_force_freq1.legend([l0, l1, l2],labels=line_labels,loc="center right", fontsize=9)
                fig_force_freq1.subplots_adjust(right=0.78)

                fig_force_freq1.savefig('freq_forces'+str(freq)+'.png')
            freq = 2
            figure_size = (8,3.55) # inches
            plt.figure(21+freq, figsize=figure_size, dpi=150) #,figsize=(6.0, 4.0), dpi=80)
            if (group==0):
                # create subplot
                fig_force_freq2, ax_force_freq2 =plt.subplots(nrows=1, ncols=2, sharey='row', squeeze=True, figsize=figure_size, dpi=150)
                fig_force_freq2.subplots_adjust(hspace=0.05)
                fig_force_freq2.subplots_adjust(wspace=0.05)
                ymin = 0.1
                ymax = 5
                ymax_pend = 1.5
                x1, x2, y1, y2 = 1.2, 2.75, 0.19, 1.75 # specify the limits for inset axis
                zoom = 2.5

                plot_num = 0
                l0 = ax_force_freq2[plot_num].plot([freq_pendulum[freq],freq_pendulum[freq]],[ymin, ymax_pend],linestyle=':',color='k')
                l1 = ax_force_freq2[plot_num].plot(w[wi_0:w_len],A_Fx[freq,wi_0:w_len],linestyle=linestyles[group],color=colors[freq],label=('forces '+hapticforces[group]))

                # set plot parameters
                ax_force_freq2[plot_num].set_xscale('log')
                ax_force_freq2[plot_num].set_yscale('log')
                # ax_force_freq2[plot_num].set_xlim((10^-1,w[w_len-1]))
                # ax_force_freq2[plot_num].set_ylim(ymin,ymax)
                ax_force_freq2[plot_num].grid(True, color=colors[4])
                ax_force_freq2[plot_num].text(.15,.93,'X-Component',horizontalalignment='left',transform=ax_force_freq2[plot_num].transAxes, fontsize=10)
                ax_force_freq2[plot_num].text(.05,.93,'A.',horizontalalignment='left',transform=ax_force_freq2[plot_num].transAxes, fontsize=10, fontweight='bold')

                plot_num = 1
                ax_force_freq2[plot_num].plot([freq_pendulum[freq],freq_pendulum[freq]],[ymin, ymax_pend],linestyle=':',color='k')
                ax_force_freq2[plot_num].plot(w[wi_0:w_len],A_Fy[freq,wi_0:w_len],linestyle=linestyles[group],color=colors[freq],label=('forces '+hapticforces[group]))

                # set plot parameters
                ax_force_freq2[plot_num].set_xscale('log')
                ax_force_freq2[plot_num].set_yscale('log')
                # ax_force_freq2[plot_num].set_xlim((10^-1,w[w_len-1]))
                # ax_force_freq2[plot_num].set_ylim(ymin,ymax)
                ax_force_freq2[plot_num].grid(True, color=colors[4])
                ax_force_freq2[plot_num].text(.15,.93,'Y-Component',horizontalalignment='left',transform=ax_force_freq2[plot_num].transAxes, fontsize=10)
                ax_force_freq2[plot_num].text(.05,.93,'B.',horizontalalignment='left',transform=ax_force_freq2[plot_num].transAxes, fontsize=10, fontweight='bold')

            if (group==1):
                plot_num = 0
                l5 = ax_force_freq2[plot_num].plot(w[wi_0:w_len],A_Fx[freq,wi_0:w_len],linestyle=linestyles[group],color=colors[freq],label=('forces '+hapticforces[group]))
                for label in (ax_force_freq2[plot_num].get_xticklabels() + ax_force_freq2[plot_num].get_yticklabels()):
                    label.set_fontsize(8)

                plot_num = 1
                ax_force_freq2[plot_num].plot(w[wi_0:w_len],A_Fy[freq,wi_0:w_len],linestyle=linestyles[group],color=colors[freq],label=('forces '+hapticforces[group]+' '+frequencylabels[freq]))
                for label in (ax_force_freq2[plot_num].get_xticklabels() + ax_force_freq2[plot_num].get_yticklabels()):
                    label.set_fontsize(8)

                # create titles
                fig_force_freq2.text(0.5, 0.91, 'Force Frequency Spectrums for 1.5Hz Ball', ha='center', fontsize=10, fontweight='bold')
                # freq.suptitle('Force Frequency Spectrums', fontsize=10, fontweight='bold')
                fig_force_freq2.text(0.5, 0.01, 'Frequency (Hz)', ha='center', fontsize=10)
                fig_force_freq2.text(0.065, 0.5, 'Normalized Frequency Amplitude', va='center', rotation='vertical', fontsize=10)
                # freq.legend(loc="upper right", fontsize=9)
                line_labels = ["Ball's resonant\nfrequency",
                    'forces '+hapticforces[0],
                    'forces '+hapticforces[1]]
                fig_force_freq2.legend([l0, l1, l2],labels=line_labels,loc="center right", fontsize=9)
                fig_force_freq2.subplots_adjust(right=0.78)

                fig_force_freq2.savefig('freq_forces'+str(freq)+'.png')
            freq = 3
            figure_size = (8,3.55) # inches
            plt.figure(21+freq, figsize=figure_size, dpi=150) #,figsize=(6.0, 4.0), dpi=80)
            if (group==0):
                # create subplot
                fig_force_freq3, ax_force_freq3 =plt.subplots(nrows=1, ncols=2, sharey='row', squeeze=True, figsize=figure_size, dpi=150)
                fig_force_freq3.subplots_adjust(hspace=0.05)
                fig_force_freq3.subplots_adjust(wspace=0.05)
                ymin = 0.1
                ymax = 5
                ymax_pend = 1.5
                x1, x2, y1, y2 = 1.2, 2.75, 0.19, 1.75 # specify the limits for inset axis
                zoom = 2.5

                plot_num = 0
                l0 = ax_force_freq3[plot_num].plot([freq_pendulum[freq],freq_pendulum[freq]],[ymin, ymax_pend],linestyle=':',color='k')
                l1 = ax_force_freq3[plot_num].plot(w[wi_0:w_len],A_Fx[freq,wi_0:w_len],linestyle=linestyles[group],color=colors[freq],label=('forces '+hapticforces[group]))

                # set plot parameters
                ax_force_freq3[plot_num].set_xscale('log')
                ax_force_freq3[plot_num].set_yscale('log')
                # ax_force_freq3[plot_num].set_xlim((10^-1,w[w_len-1]))
                # ax_force_freq3[plot_num].set_ylim(ymin,ymax)
                ax_force_freq3[plot_num].grid(True, color=colors[4])
                ax_force_freq3[plot_num].text(.15,.93,'X-Component',horizontalalignment='left',transform=ax_force_freq3[plot_num].transAxes, fontsize=10)
                ax_force_freq3[plot_num].text(.05,.93,'A.',horizontalalignment='left',transform=ax_force_freq3[plot_num].transAxes, fontsize=10, fontweight='bold')

                plot_num = 1
                ax_force_freq3[plot_num].plot([freq_pendulum[freq],freq_pendulum[freq]],[ymin, ymax_pend],linestyle=':',color='k')
                ax_force_freq3[plot_num].plot(w[wi_0:w_len],A_Fy[freq,wi_0:w_len],linestyle=linestyles[group],color=colors[freq],label=('forces '+hapticforces[group]))

                # set plot parameters
                ax_force_freq3[plot_num].set_xscale('log')
                ax_force_freq3[plot_num].set_yscale('log')
                # ax_force_freq3[plot_num].set_xlim((10^-1,w[w_len-1]))
                # ax_force_freq3[plot_num].set_ylim(ymin,ymax)
                ax_force_freq3[plot_num].grid(True, color=colors[4])
                ax_force_freq3[plot_num].text(.15,.93,'Y-Component',horizontalalignment='left',transform=ax_force_freq3[plot_num].transAxes, fontsize=10)
                ax_force_freq3[plot_num].text(.05,.93,'B.',horizontalalignment='left',transform=ax_force_freq3[plot_num].transAxes, fontsize=10, fontweight='bold')

            if (group==1):
                plot_num = 0
                l5 = ax_force_freq3[plot_num].plot(w[wi_0:w_len],A_Fx[freq,wi_0:w_len],linestyle=linestyles[group],color=colors[freq],label=('forces '+hapticforces[group]))
                for label in (ax_force_freq3[plot_num].get_xticklabels() + ax_force_freq3[plot_num].get_yticklabels()):
                    label.set_fontsize(8)

                plot_num = 1
                ax_force_freq3[plot_num].plot(w[wi_0:w_len],A_Fy[freq,wi_0:w_len],linestyle=linestyles[group],color=colors[freq],label=('forces '+hapticforces[group]+' '+frequencylabels[freq]))
                for label in (ax_force_freq3[plot_num].get_xticklabels() + ax_force_freq3[plot_num].get_yticklabels()):
                    label.set_fontsize(8)

                # create titles
                fig_force_freq3.text(0.5, 0.91, 'Force Frequency Spectrums for 2Hz Ball', ha='center', fontsize=10, fontweight='bold')
                # freq.suptitle('Force Frequency Spectrums', fontsize=10, fontweight='bold')
                fig_force_freq3.text(0.5, 0.01, 'Frequency (Hz)', ha='center', fontsize=10)
                fig_force_freq3.text(0.065, 0.5, 'Normalized Frequency Amplitude', va='center', rotation='vertical', fontsize=10)
                # freq.legend(loc="upper right", fontsize=9)
                line_labels = ["Ball's resonant\nfrequency",
                    'forces '+hapticforces[0],
                    'forces '+hapticforces[1]]
                fig_force_freq3.legend([l0, l1, l2],labels=line_labels,loc="center right", fontsize=9)
                fig_force_freq3.subplots_adjust(right=0.78)

                fig_force_freq3.savefig('freq_forces'+str(freq)+'.png')


plt.show()
