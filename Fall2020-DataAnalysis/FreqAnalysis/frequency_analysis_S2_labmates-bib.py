import numpy as np
import scipy as sp
from numpy import fft, genfromtxt
from scipy import fft, arange, stats
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
# from make_boxplot import make_boxplot, add_stats
from plot_utils import *
from plot_freq_spectrum import xy_spectrum,mag_spectrum
import csv
from perform_transform import calculate_amplitude,normalize_spectrum

################################################################################
# This program analyses the frequency spectrum for each trial/participant.
# It saves and shows figures for the BioRob paper and saves the values of the peaks
# around the pendulum's resonant frequency into csv files
# "freq-metrics.csv", "freq-metrics-controls.csv", "freq-metrics-stroke.csv"
# for statistical analyses in R.
# This code expects the original data to be formatted and named properly
################################################################################

# Edit these variables before running
make_plots = 1 # 0-do not make plots 1-make plots
make_plot_each_sub = 0 # 0-do not make plots 1-make plots
# only_noforces = 0
haptic_forces_added = 0
freq_range_boxplot = 1
window = .2


number_of_subjects = 3
min_sub = 1
DIR = "Z:Fall2020Data-round2/" #set directory where data is mounted Ola- "/media/ola/Elements/R01prelim" Milli -"Z:"
DT = 0.05
Fs = 1/DT
freq_pendulum = [.5,1,1.5,2.5]

# Label factors as strings for ezANOVA analysis
# forceconditions = ['baseline','forces no ball','no forces','forces']
forceconditions = ['Baseline Motion','Still Ball With Haptic Forces','Moving Ball Without Haptic Forces','Moving Ball With Haptic Forces']
group_label = ['F0_B0','F1_B0','F0_B1','F1_B1']
frequencylabels = ['0.5Hz','1Hz','1.5Hz','2.5Hz']
subjects = ['Subject1','Subject2','Subject3','Subject4','Subject5','Subject6','Subject7','Subject8','Subject9','Subject10','Subject11','Subject12','Subject13','Subject14','Subject15','Subject16','Subject17','Subject18']
# colors = ['#601A4A','#EE442F','#63ACBE','#006400'] # for plots
# linestyles = ['--','-']
markerstyles = ['o','D']
# radius_options = [ 0.995, 0.249, 0.111, 0.062, 0.04 ];

# Create vector of frequencies of interest
nyquist_freq = int(np.floor(Fs/2))
# nyquist_freq = 40 # override the nyquist frequency (there is little info of value at higher frequencies)
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

# create matrices for saving energy
# if freq_range_boxplot:
#     window = .2 # for calculating the energy around different frequencies
#     energy_mat_fx = np.zeros((4,4,number_of_subjects))
#     energy_mat_fy = np.zeros((4,4,number_of_subjects))
#     energy_num = np.zeros((4,4,number_of_subjects)) # saves the number of times trials for each condition

if make_plots:
    xy_xlist = []
    xy_ylist = []
    mag_list = []
# if plotting each subject, iterate through each subject plot
# num_iters = 1
# if make_plot_each_sub or freq_range_boxplot:
#     num_iters = number_of_subjects
# for sub_plot in range(1,1+num_iters):

# Iterate through all the files starting with with and without haptic forces
compare_groups = np.zeros((12,number_of_subjects))
for group in range(len(forceconditions)):
    # if freq_range_boxplot:
    #     if group==0 or group==1 or group==2:
    #         continue
    force_str = group_label[group]

    num_freq = np.zeros(4)
    A_Fmag = np.zeros((4,w_len))
    A_Fx = np.zeros((4,w_len))
    A_Fy = np.zeros((4,w_len))

    if freq_range_boxplot:
        energy_mat_fx = np.zeros((4,4,number_of_subjects))
        energy_mat_fy = np.zeros((4,4,number_of_subjects))
        energy_num = np.zeros((4,4,number_of_subjects)) # saves the number of times trials for each condition

    if haptic_forces_added:
        print('Initialize haptic force trials')
        A_Fx_forceadd = np.zeros((4,w_len))
        A_Fy_forceadd = np.zeros((4,w_len))
        if group==0:
            print('Skip no haptic force trials')
            continue

    for subject_num in range(min_sub,1+number_of_subjects): #iterate though subjects

        if (make_plot_each_sub==0) or (make_plot_each_sub==1 and sub_plot==subject_num):
            # Sets the correct subfile
            subname = "S0" + str(subject_num)
            subfile = "/OutputData_S" + str(subject_num)

            for freq in range(0,4):
                for k in range(1,80): #iterate through trials
                    try:

                        # Open the trial files
                        trialfile = DIR+subname+subfile+"_Trial"+str(k)+"_Freq"+str(freq)+"_SL0_"+force_str+".csv"
                        data = genfromtxt(trialfile,delimiter=',',dtype=float)
                        data = np.delete(data,0,0) # Deletes the first column of column names
                        print(trialfile)


                        # Get amplitudes for force and position signals
                        # Combine the x and y direction forces
                        y = np.sqrt(np.square(data[:,8])+np.square(data[:,9]))
                        # y = np.square(data[:,8])+np.square(data[:,9])
                        A_Fmag_i,frq = calculate_amplitude(w,y,Fs,True,False)

                        # Use signal directly 1-x position 2- y position 8- x force 9- y force 12-x ball force 13-y ball force
                        if haptic_forces_added: # add the haptic forces
                            y_Fx = data[:,8]+data[:,12]
                            A_Fx_i,frq = calculate_amplitude(w,y_Fx,Fs,False,True)
                            y_Fy = data[:,9]+data[:,13]
                            A_Fy_i,frq = calculate_amplitude(w,y_Fy,Fs,False,True)
                            A_Fx_forceadd[freq,:] += A_Fx_i
                            A_Fy_forceadd[freq,:] += A_Fy_i
                        y_Fx = data[:,8]
                        A_Fx_i,frq = calculate_amplitude(w,y_Fx,Fs,True,False)
                        y_Fy = data[:,9]
                        A_Fy_i,frq = calculate_amplitude(w,y_Fy,Fs,True,False)
                        A_Fmag[freq,:] += A_Fmag_i
                        A_Fx[freq,:] += A_Fx_i
                        A_Fy[freq,:] += A_Fy_i
                        num_freq[freq] += 1


                        if freq_range_boxplot:
                            dw = w[1]-w[0]
                            for freq2 in range(0,4):
                                Ax_norm = normalize_spectrum(A_Fx[freq,:],dw)
                                Ay_norm = normalize_spectrum(A_Fy[freq,:],dw)

                                freq_list = []
                                for w_i in range(0,len(w)):
                                    if (w[w_i] < freq_pendulum[freq2]+window) and (w[w_i] > freq_pendulum[freq2]-window):
                                        freq_list.append(Ax_norm[w_i])
                                energy_mat_fx[freq,freq2,subject_num-1] += np.sum(np.square(freq_list))*dw
                                # print(energy_mat_fx[freq,freq2,subject_num-1])
                                freq_list = []
                                for w_i in range(0,len(w)):
                                    if (w[w_i] < freq_pendulum[freq2]+window) and (w[w_i] > freq_pendulum[freq2]-window):
                                        freq_list.append(Ay_norm[w_i])
                                energy_mat_fy[freq,freq2,subject_num-1] += np.sum(np.square(freq_list))*dw

                                energy_num[freq,freq2,subject_num-1] += 1
                    except:
                        pass

    if make_plots==1:
        # Take average of each signal

        for freq in range(0,4):
            if group == 0:
                if freq > 0:
                    A_Fmag[0,:] += A_Fmag[freq,:]
                    A_Fx[0,:] += A_Fx[freq,:]
                    A_Fy[0,:] += A_Fy[freq,:]
                    num_freq[0] += num_freq[freq]
            else:
                if num_freq[freq]==0:
                    print('No data was added for group', group, 'frequency', freq)
                else:
                    dw = w[1]-w[0]
                    A_Fmag[freq,:] /= num_freq[freq]
                    A_Fx[freq,:] /= num_freq[freq]
                    A_Fy[freq,:] /= num_freq[freq]
                    A_Fmag[freq,:] = normalize_spectrum(A_Fmag[freq,:],dw)
                    A_Fx[freq,:] = normalize_spectrum(A_Fx[freq,:],dw)
                    A_Fy[freq,:] = normalize_spectrum(A_Fy[freq,:],dw)
                    mag_list.append(A_Fmag[freq,:])
                    xy_xlist.append(A_Fx[freq,:])
                    xy_ylist.append(A_Fy[freq,:])

                    # if freq_range_boxplot and group==0:
                    #     dw = w[1]-w[0]
                    #     for freq2 in range(0,4):
                    #         freq_list = []
                    #         for w_i in range(0,len(w)):
                    #             if (w[w_i] < freq_pendulum[freq2]+window) and (w[w_i] > freq_pendulum[freq2]-window):
                    #                 freq_list.append(A_Fx[freq,w_i])
                    #         energy_mat_fx[freq,freq2,sub_plot-1] += np.sum(np.square(freq_list))*dw
                    #
                    #         freq_list = []
                    #         for w_i in range(0,len(w)):
                    #             if (w[w_i] < freq_pendulum[freq2]+window) and (w[w_i] > freq_pendulum[freq2]-window):
                    #                 freq_list.append(A_Fy[freq,w_i])
                    #         energy_mat_fy[freq,freq2,sub_plot-1] += np.sum(np.square(freq_list))*dw
        if group == 0:
            freq = 0
            dw = w[1]-w[0]
            A_Fmag[freq,:] /= num_freq[freq]
            A_Fx[freq,:] /= num_freq[freq]
            A_Fy[freq,:] /= num_freq[freq]
            A_Fmag[freq,:] = normalize_spectrum(A_Fmag[freq,:],dw)
            A_Fx[freq,:] = normalize_spectrum(A_Fx[freq,:],dw)
            A_Fy[freq,:] = normalize_spectrum(A_Fy[freq,:],dw)

        for freq in range(0,4):
            if group == 0:
                mag_list.append(A_Fmag[0,:])
                xy_xlist.append(A_Fx[0,:])
                xy_ylist.append(A_Fy[0,:])
            if haptic_forces_added:
                # A_Fmag_forceadd[freq,:] /= num_freq[freq]
                A_Fx_forceadd[freq,:] /= num_freq[freq]
                A_Fy_forceadd[freq,:] /= num_freq[freq]
                # A_Fmag_forceadd[freq,:] = normalize_spectrum(A_Fmag_forceadd[freq,:],dw)
                A_Fx_forceadd[freq,:] = normalize_spectrum(A_Fx_forceadd[freq,:],dw)
                A_Fy_forceadd[freq,:] = normalize_spectrum(A_Fy_forceadd[freq,:],dw)

                xy_xlist.append(A_Fx_forceadd[freq,:])
                xy_ylist.append(A_Fy_forceadd[freq,:])


        # Plot the frequency spectrum for this group
        title = 'Frequency Spectrum Plot: ' + forceconditions[group]
        xlabel = 'Frequency (Hz)'
        ylabel = 'Normalized Frequency Amplitude'
        colors = ['#601A4A','#EE442F','#63ACBE','#006400']
        linestyles = ['-','-','-','-']
        ymin = 0.1
        ymax = 2
        ymax_pend = 1.5
        if group==0:
            legend = ["Ball's resonant\nfrequency"]
        else:
            legend = ["Ball's resonant\nfrequency",
                frequencylabels[0],
                frequencylabels[1],
                frequencylabels[2],
                frequencylabels[3]]
        n = len(xy_xlist)
        [fig,ax] = xy_spectrum(w,xy_xlist[n-4:n],xy_ylist[n-4:n],freq_pendulum,title,xlabel,ylabel,legend,linestyles,colors,ymin,ymax,ymax_pend)
        fig.savefig('Plots/'+'xy_freq_'+group_label[group]+'.png')

        # Plot boxplots for this group
        if group != 0:
            energy_mat_fx = np.divide(energy_mat_fx,energy_num)
            energy_mat_fy = np.divide(energy_mat_fy,energy_num)
            data = []
            for freq2 in range(0,4):
                for freq in range(0,4):
                    energy = (energy_mat_fx[freq,freq2,:] + energy_mat_fy[freq,freq2,:])/2
                    data.append(np.array(energy))

            # create lists for colors, labels, ect.
            n = len(data)
            labels = []
            box_alpha = []
            box_colors = []
            for i in range(n):
                labels.append(freq_pendulum[i % 4])
                box_alpha.append(1)
                if (i // 4) == (i % 4): # if the ball and range match
                    box_colors.append('#ff7e0d')
                    compare_groups[int((3*(i//4))+group-1),:] = data[i]
                else:
                    box_colors.append('#d6caed')


            figure_size = (6,3.55) # sets the size of the figure in inches
            xlabel = ''
            ylabel = 'Fraction of Total Energy'

            title = 'Energy Content of Movement at Different Frequencies:\n' + forceconditions[group]
            [fig,ax]=make_boxplot(data,title,xlabel,ylabel,labels,box_colors,box_alpha,figure_size)

            # Add rectangles for background
            [ymin,ymax]=ax.get_ylim()
            grey_colors = ['#e3e3e3','#ffffff','#e3e3e3','#ffffff']
            x1 = []
            x2 = []
            for i in range(4):
                x = i*4+0.5
                ax.add_patch(Rectangle((x,ymin), 4, ymax-ymin, facecolor=grey_colors[i], alpha=1,zorder=1))
                x1.append(x)
                if i!=0:
                    x2.append(x)
            x2.append(x+4)

            # Text Labels
            y = ymin - ((ymax-ymin)/8)
            text_buffer = ((ymax-ymin)/20)
            name = ['0.5+/-0.2','1+/-0.2','1.5+/-0.2','2.5+/-0.2']
            add_labels(ax,x1,x2,y,name,text_buffer)
            y = ymin - ((ymax-ymin)/6)
            ax.text(0.5, y,'Range(Hz):',horizontalalignment='right', fontname="Arial", fontsize=10, fontweight='bold')
            y = ymin - ((ymax-ymin)/12)
            ax.text(0.5, y,'Ball(Hz):',horizontalalignment='right', fontname="Arial", fontsize=10, fontweight='bold')
            # add_labels(ax,x1,x2,y,name,text_buffer)
            fig.savefig('Plots/'+'energy_'+group_label[group]+'.png')


if make_plots:

    ####################################################################
    ####################################################################
    ######   Compare the different haptic force conditions  #####
    ####################################################################
    ####################################################################

    figure_size = (7,3) # sets the size of the figure in inches
    xlabel = ''
    ylabel = 'Fraction of Total Energy'

    title = 'Energy Content of Movement at Resonance'

    n = compare_groups.shape[0]
    data = []
    labels = []
    box_alpha = []
    box_colors = []
    for i in range(n):
        data.append(np.array(compare_groups[i,:]))
        labels.append('')
        box_alpha.append(1)
        box_colors.append(colors[i % 3])

    [fig,ax]=make_boxplot(data,title,xlabel,ylabel,labels,box_colors,box_alpha,figure_size)

    # Add rectangles for background
    [ymin,ymax]=ax.get_ylim()
    grey_colors = ['#e3e3e3','#ffffff','#e3e3e3','#ffffff']
    x1 = []
    x2 = []
    for i in range(4):
        x = i*3+0.5
        ax.add_patch(Rectangle((x,ymin), 3, ymax-ymin, facecolor=grey_colors[i], alpha=1,zorder=1))
        x1.append(x)
        if i!=0:
            x2.append(x)
    x2.append(x+3)

    # Text Labels
    y = ymin - ((ymax-ymin)/15)
    text_buffer = ((ymax-ymin)/15)
    name =frequencylabels
    add_labels(ax,x1,x2,y,name,text_buffer)
    y = ymin - ((ymax-ymin)/7.5)
    ax.text(0.5, y,'Ball:',horizontalalignment='right', fontname="Arial", fontsize=10, fontweight='bold')

    x = 5
    l1 = ax.add_patch(Rectangle((x,ymin-100), 3, ymax-ymin, facecolor=colors[0], alpha=1,zorder=1))
    l2 = ax.add_patch(Rectangle((x,ymin-100), 3, ymax-ymin, facecolor=colors[1], alpha=1,zorder=1))
    l3 = ax.add_patch(Rectangle((x,ymin-100), 3, ymax-ymin, facecolor=colors[2], alpha=1,zorder=1))
    ax.set_ylim(ymin,ymax)
    legend_cond = ['Still Ball With\nHaptic Forces','Moving Ball Without\nHaptic Forces','Moving Ball With\nHaptic Forces']

    fig.subplots_adjust(right=0.72)
    L = fig.legend([l1,l2,l3], legend_cond,ncol=1, fontsize=10,loc='center right')

    fig.savefig('Plots/'+'energy_comparegroups.png')

    ####################################################################
    ####################################################################
    ######   Frequency spectrum plot for each frequency  #####
    ####################################################################
    ####################################################################
    title = 'Frequency Spectrum Plot'
    xlabel = 'Frequency (Hz)'
    ylabel = 'Normalized Frequency Amplitude'
    colors = ['#601A4A','#EE442F','#63ACBE','#006400']
    linestyles = ['-','-','-','-']
    ymin = 0.1
    ymax = 2
    ymax_pend = 1.5
    for freq in range(0,4):
        title = 'Frequency Spectrum Plot for ' + frequencylabels[freq]
        savename = 'xy_freq_'+frequencylabels[freq]+'.png'
        legend = ["Ball's resonant\nfrequency",
                    'Baseline Motion',
                    'Still Ball With\nHaptic Forces',
                    'Moving Ball Without\nHaptic Forces',
                    'Moving Ball With\nHaptic Forces']
        # # colors_fq = [colors[freq],colors[freq+4]]
        # colors_fq = [colors[0],colors[1],colors[2],colors[3]]
        # # linestyles_fq = [linestyles[freq],linestyles[freq+4]]
        # linestyles_fq = [linestyles[4],linestyles[5],linestyles[6],linestyles[7]]
        xdata = []
        ydata = []
        xdata.append(xy_xlist[freq])
        xdata.append(xy_xlist[freq+4])
        xdata.append(xy_xlist[freq+8])
        xdata.append(xy_xlist[freq+12])
        ydata.append(xy_ylist[freq])
        ydata.append(xy_ylist[freq+4])
        ydata.append(xy_ylist[freq+8])
        ydata.append(xy_ylist[freq+12])
        [fig,ax] = xy_spectrum(w,xdata,ydata,[freq_pendulum[freq]],title,xlabel,ylabel,legend,linestyles,colors,ymin,ymax,ymax_pend)
        fig.savefig('Plots/'+savename)

plt.show()
