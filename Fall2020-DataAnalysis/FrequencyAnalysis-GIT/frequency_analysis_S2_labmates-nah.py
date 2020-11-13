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
from plot_freq_spectrum import xy_spectrum,mag_spectrum

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
number_of_subjects = 5 # 7
DIR = "/media/ola/Data/Research/DowntownData/Fall2020LabData/"
DT = 0.05
Fs = 1/DT
freq_pendulum = [.5,1,1.5,2]
side_arm = ["R","R","L","L","R","L"]
split_signal = 1

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
frequencylabels = ['R','L']
subjects = ['Subject1','Subject2','Subject3','Subject4','Subject5','Subject6','Subject7','Subject8','Subject9','Subject10','Subject11','Subject12','Subject13','Subject14','Subject15','Subject16','Subject17','Subject18']
#colors = ['#601A4A','#EE442F','#63ACBE','#006400',"#E0E0E0"] # for plots
#linestyles = ['--','-']
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

if make_plots:
    xy_xlist = []
    xy_ylist = []
    mag_list = []

for sub_plot in range(2,1+num_iters):
    # if make_plot_each_sub:
    #     # For plotting individual subject figures
    #     plt.figure(plot_num)
    #     plt.title('Subject '+str(plot_num)+' Frequency Spectrum')

    # Iterate through all the files starting with with and without haptic forces
    for group in range(0,1):

        # The average of each amplitude signal is taken. num_0percent-num_50percent records the number of signals added together
        num_freq = np.zeros(2)
        A_Fmag = np.zeros((2,w_len))
        A_Fx = np.zeros((2,w_len))
        A_Fy = np.zeros((2,w_len))
        if split_signal:
            A_Fmag_nails = np.zeros((2,w_len))
            A_Fx_nails = np.zeros((2,w_len))
            A_Fy_nails = np.zeros((2,w_len))
            A_Fmag_other = np.zeros((2,w_len))
            A_Fx_other = np.zeros((2,w_len))
            A_Fy_other = np.zeros((2,w_len))
        for subject_num in range(2,1+number_of_subjects): #iterate though subjects

            if (make_plot_each_sub==0) or (make_plot_each_sub==1 and sub_plot==subject_num):
                # Sets the correct subfile
                subname = "S0" + str(subject_num)
                subfile = "/NailHammerData_S" + str(subject_num)

                for side in range(0,2): # iterate through R/L
                    for k in range(1,13): #iterate through trials
                        try:
                            # Open the trial files
                            trialfile = DIR+subname+subfile+"_Trial"+str(k)+"_Task"+str(side)+"_SL0.csv"
                            # print(trialfile)
                            data = genfromtxt(trialfile,delimiter=',',dtype=float)
                            data = np.delete(data,0,0) # Deletes the first column of column names
                            print(trialfile)

                            if ((side_arm[subject_num]=="R" and side==0) or (side_arm[subject_num]=="L" and side==1)):
                                freq = 0
                            else:
                                freq = 1

                            # Get amplitudes for force and position signals
                            # Combine the x and y direction forces
                            y = np.sqrt(np.square(data[:,9])+np.square(data[:,10]))
                            A_Fmag_i,frq = calculate_amplitude(w,y,Fs,False,True)
                            # print(A_Fmag_i)
                            # isofhaihs
                            # Use signal directly 1-x position 2- y position 8- x force 9- y force
                            y_Fx = data[:,9]
                            A_Fx_i,frq = calculate_amplitude(w,y_Fx,Fs,False,True)
                            y_Fy = data[:,10]
                            A_Fy_i,frq = calculate_amplitude(w,y_Fy,Fs,False,True)

                            if make_plots==1:
                                # For combining participants for plots
                                A_Fmag[freq,:] += A_Fmag_i
                                A_Fx[freq,:] += A_Fx_i
                                A_Fy[freq,:] += A_Fy_i
                                num_freq[freq] += 1

                            if split_signal:
                                y_nails = []; y_Fx_nails = []; y_Fy_nails = []
                                y_other = []; y_Fx_other = []; y_Fy_other = []
                                for row in range(0,data.shape[0]):
                                    if data[row,12] % 5>0:
                                        y_nails.append(np.sqrt(np.square(data[row,9])+np.square(data[row,10])))
                                        y_Fx_nails.append(data[row,9])
                                        y_Fy_nails.append(data[row,10])
                                    elif data[row,12] % 5==0:
                                        y_other.append(np.sqrt(np.square(data[row,9])+np.square(data[row,10])))
                                        y_Fx_other.append(data[row,9])
                                        y_Fy_other.append(data[row,10])
                                A_Fmag_nails_i,frq = calculate_amplitude(w,y_nails,Fs,False,True)
                                A_Fx_nails_i,frq = calculate_amplitude(w,y_Fx_nails,Fs,False,True)
                                A_Fy_nails_i,frq = calculate_amplitude(w,y_Fy_nails,Fs,False,True)

                                A_Fmag_other_i,frq = calculate_amplitude(w,y_other,Fs,False,True)
                                A_Fx_other_i,frq = calculate_amplitude(w,y_Fx_other,Fs,False,True)
                                A_Fy_other_i,frq = calculate_amplitude(w,y_Fy_other,Fs,False,True)

                                if np.min(A_Fmag_nails_i)>0:
                                    A_Fmag_nails[freq,:] += A_Fmag_nails_i
                                    A_Fx_nails[freq,:] += A_Fx_nails_i
                                    A_Fy_nails[freq,:] += A_Fy_nails_i

                                if np.min(A_Fmag_other_i)>0:
                                    A_Fmag_other[freq,:] += A_Fmag_other_i
                                    A_Fx_other[freq,:] += A_Fx_other_i
                                    A_Fy_other[freq,:] += A_Fy_other_i

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
            for freq in range(0,2):
                if num_freq[freq]!=0:
                    dw = w[1]-w[0]
                    print("dw ", dw)
                    if split_signal:
                        A_Fmag_nails[freq,:] /= num_freq[freq]
                        A_Fx_nails[freq,:] /= num_freq[freq]
                        A_Fy_nails[freq,:] /= num_freq[freq]
                        E_signal = np.sum(np.square(A_Fmag_nails[freq,:]))*dw
                        A_Fmag_nails[freq,:] *= np.sqrt(1/E_signal)
                        E_signal = np.sum(np.square(A_Fx_nails[freq,:]))*dw
                        A_Fx_nails[freq,:] *= np.sqrt(1/E_signal)
                        E_signal = np.sum(np.square(A_Fy_nails[freq,:]))*dw
                        A_Fy_nails[freq,:] *= np.sqrt(1/E_signal)

                        A_Fmag_other[freq,:] /= num_freq[freq]
                        A_Fx_other[freq,:] /= num_freq[freq]
                        A_Fy_other[freq,:] /= num_freq[freq]
                        E_signal = np.sum(np.square(A_Fmag_other[freq,:]))*dw
                        A_Fmag_other[freq,:] *= np.sqrt(1/E_signal)
                        E_signal = np.sum(np.square(A_Fx_other[freq,:]))*dw
                        A_Fx_other[freq,:] *= np.sqrt(1/E_signal)
                        E_signal = np.sum(np.square(A_Fy_other[freq,:]))*dw
                        A_Fy_other[freq,:] *= np.sqrt(1/E_signal)
                    else:
                        A_Fmag[freq,:] /= num_freq[freq]
                        A_Fx[freq,:] /= num_freq[freq]
                        A_Fy[freq,:] /= num_freq[freq]

            if split_signal:# and group==1:
                # print(np.shape(A_Fmag_other[0,:]))
                # print(np.shape(A_Fmag_other[1,:]))
                # print(A_Fmag_other[1,:])
                for freq in range(0,2):
                    mag_list.append(A_Fmag_other[freq,:])
                    xy_xlist.append(A_Fx_other[freq,:])
                    xy_ylist.append(A_Fy_other[freq,:])
                    mag_list.append(A_Fmag_nails[freq,:])
                    xy_xlist.append(A_Fx_nails[freq,:])
                    xy_ylist.append(A_Fy_nails[freq,:])
                # print(np.shape(mag_list))
                # print(np.shape(xy_xlist))
                # print(np.shape(xy_ylist))


if make_plots:
    ####################################################################
    ####################################################################
    ###############   Frequency spectrum plot (x/y)    #################
    ####################################################################
    ####################################################################
    title = 'Frequency Spectrum Plot'
    savename = 'xy_freq.png'
    xlabel = 'Frequency (Hz)'
    ylabel = 'Normalized Frequency Amplitude'
    colors = ['#601A4A','#EE442F','#63ACBE','#006400','#601A4A','#EE442F','#63ACBE','#006400']
    linestyles = ['--','--','--','--','-','-','-','-']
    if split_signal:
        legend = ["",
            'Non-nail movement\n(same side)',
            'Non-nail movement\n(cross side)',
            'Nail movement\n(same side)',
            'Nail movement\n(cross side)']
        savename = 'xy_freq_split.png'
    # else:
    #     legend = ["Ball's resonant\nfrequency",
    #         'forces '+hapticforces[0]+' '+frequencylabels[0],
    #         'forces '+hapticforces[0]+' '+frequencylabels[1],
    #         'forces '+hapticforces[0]+' '+frequencylabels[2],
    #         'forces '+hapticforces[0]+' '+frequencylabels[3],
    #         'forces '+hapticforces[1]+' '+frequencylabels[0],
    #         'forces '+hapticforces[1]+' '+frequencylabels[1],
    #         'forces '+hapticforces[1]+' '+frequencylabels[2],
    #         'forces '+hapticforces[1]+' '+frequencylabels[3]]
        # savename = 'xy_freq.png'

    ymin = 0.1
    ymax = 2
    ymax_pend = 1.5
    [fig,ax] = xy_spectrum(w,xy_xlist,xy_ylist,freq_pendulum,title,xlabel,ylabel,legend,linestyles,colors,ymin,ymax,ymax_pend)
    fig.savefig(savename)

    ####################################################################
    ####################################################################
    ##############   Frequency magnitude spectrum plot  ################
    ####################################################################
    ####################################################################
    title = 'Frequency Magnitude Spectrum Plot'
    if split_signal:
        savename = 'mag_freq_split.png'
    else:
        savename = 'mag_freq.png'
    ymin = 0.01
    [fig,ax] = mag_spectrum(w,mag_list,freq_pendulum*2,title,xlabel,ylabel,legend,linestyles,colors,ymin,ymax,ymax_pend)
    fig.savefig(savename)

    ####################################################################
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
        l1 = ax_force[plot_num].plot(w[wi_0:w_len],A_Fx[0,wi_0:w_len],linestyle=linestyles[group],color=colors[0],label=frequencylabels[0])
        l2 = ax_force[plot_num].plot(w[wi_0:w_len],A_Fx[1,wi_0:w_len],linestyle=linestyles[group],color=colors[1],label=frequencylabels[1])

        # set plot parameters
        ax_force[plot_num].set_xscale('log')
        ax_force[plot_num].set_yscale('log')
        #ax_force[plot_num].set_xlim((10^-1,w[w_len-1]))
        #ax_force[plot_num].set_ylim(ymin,ymax)
        ax_force[plot_num].grid(True, color=colors[4])
        ax_force[plot_num].text(.15,.93,'X-Component',horizontalalignment='left',transform=ax_force[plot_num].transAxes, fontsize=10)
        ax_force[plot_num].text(.05,.93,'A.',horizontalalignment='left',transform=ax_force[plot_num].transAxes, fontsize=10, fontweight='bold')

        plot_num = 1
        ax_force[plot_num].plot([0.5, 0.5],[ymin, ymax_pend],linestyle=':',color='k')
        for freq in range(0,2):
            ax_force[plot_num].plot(w[wi_0:w_len],A_Fy[freq,wi_0:w_len],linestyle=linestyles[group],color=colors[freq],label=frequencylabels[freq])

        # set plot parameters
        ax_force[plot_num].set_xscale('log')
        ax_force[plot_num].set_yscale('log')
        # ax_force[plot_num].set_xlim((10^-1,w[w_len-1]))
        # ax_force[plot_num].set_ylim(ymin,ymax)
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


plt.show()
