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

# Edit these variables before running
save_values = 0 # 0-do not save values for statistical tests 1-do save values
make_plots = 1 # 0-do not make plots 1-make plots
make_plot_each_sub = 1 # 0-do not make plots 1-make plots
haptic_forces_added = 0 # Include adding the haptic forces as an experimental condition
filter = 1
window = .2

number_of_subjects = 7
min_sub = 1
DIR = "Z:Fall2020Data-round2/" #set directory where data is mounted Ola- "/media/ola/Elements/R01prelim" Milli -"Z:"
DT = 0.05
Fs = 1/DT
freq_pendulum = [.5,1,1.5,2.5]

# Label factors as strings for ezANOVA analysis
forceconditions = ['Baseline Motion','Still Ball With Haptic Forces','Moving Ball Without Haptic Forces','Moving Ball With Haptic Forces']
group_label_read = ['F0_B0','F1_B0','F0_B1','F1_B1']
group_label = ['F0_B0','F1_B0','F0_B1','F1_B1']
if haptic_forces_added:
    forceconditions = ['Baseline Motion','Still Ball With Haptic Forces','Moving Ball Without Haptic Forces','Haptic Forces Added to Signal','Moving Ball With Haptic Forces']
    group_label_read = ['F0_B0','F1_B0','F0_B1','F1_B1','F1_B1']
    group_label = ['F0_B0','F1_B0','F0_B1','F1_B1_add','F1_B1']
frequencylabels = ['0.5Hz','1Hz','1.5Hz','2.5Hz']
subjects = ['Subject1','Subject2','Subject3','Subject4','Subject5','Subject6','Subject7','Subject8','Subject9','Subject10','Subject11','Subject12','Subject13','Subject14','Subject15','Subject16','Subject17','Subject18']
markerstyles = ['o','D']

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

# Store data for statistical analysis
if save_values==1:
    file_metrics = "freq-metrics.csv"
    columns = ["Subject","Force","BallFreq","Energy0.5","Energy1.0","Energy1.5","Energy2.5","Resonance"]
    with open(file_metrics,'w') as csvfile:
        testwriter = csv.writer(csvfile,delimiter=',')
        testwriter.writerow(columns)

if make_plots:
    xy_xlist = []
    xy_ylist = []
    mag_list = []

# if plotting each subject, iterate through each subject plot
num_sub_plots = 1
if make_plot_each_sub:
    num_sub_plots = number_of_subjects
for sub_plot in range(1,1+num_sub_plots):

    # Iterate through all the files starting with with and without haptic forces
    compare_groups1 = np.zeros((2,4,number_of_subjects))
    compare_groups2 = np.zeros((2,4,number_of_subjects))
    compare_groups3 = np.zeros((len(forceconditions),4,number_of_subjects))
    for group in range(0,len(forceconditions)):
        num_freq = np.zeros(4)
        A_Fmag = np.zeros((4,w_len))
        A_Fx = np.zeros((4,w_len))
        A_Fy = np.zeros((4,w_len))

        energy_mat = np.zeros((4,4,number_of_subjects))
        energy_num = np.zeros((4,4,number_of_subjects)) # saves the number of times trials for each condition
        for subject_num in range(min_sub,1+number_of_subjects): #iterate though subjects

            if (make_plot_each_sub==0) or (make_plot_each_sub==1 and sub_plot==subject_num):

                for freq in range(0,4):
                    pathName = DIR + "S0" + str(subject_num) + "/"
                    Files = []
                    fileNames = os.listdir(pathName)
                    for fileName in fileNames:
                        if fileName.startswith("OutputData"):
                            if "_Freq"+str(freq) in fileName:
                                if group_label_read[group] in fileName:
                                    Files.append(fileName)

                    for i in Files:

                        print(i)
                        file = open(os.path.join(pathName, i), "rU")
                        data = genfromtxt(file,delimiter=',',dtype=float)
                        data = np.delete(data,0,0) # Deletes the first column of column names

                        str_split = i.split("_")
                        trial_num = int(str_split[2].split("Trial")[1])
                        if subject_num==6 and trial_num==72:
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
                        for freq2 in range(0,4):
                            freq_list = []
                            for w_i in range(0,len(w)):
                                if (w[w_i] < freq_pendulum[freq2]+window) and (w[w_i] > freq_pendulum[freq2]-window):
                                    # freq_list.append(Ax_norm[w_i])
                                    # freq_list.append(Ay_norm[w_i])
                                    freq_list.append(Amag_norm[w_i])

                            energy_mat[freq,freq2,subject_num-1] += np.sum(np.square(freq_list))*dw
                            energy_num[freq,freq2,subject_num-1] += 1

                        if save_values:
                            if subject_num==6 and (trial_num==70 or trial_num==71):
                                if trial_num==70:
                                    metrics_save = [energy_mat[freq,0,subject_num-1],energy_mat[freq,1,subject_num-1],
                                                energy_mat[freq,2,subject_num-1],energy_mat[freq,3,subject_num-1]]
                                else:
                                    metrics_save =  [energy_mat[freq,0,subject_num-1]+metrics_save[0],
                                                    energy_mat[freq,1,subject_num-1]+metrics_save[1],
                                                    energy_mat[freq,2,subject_num-1]+metrics_save[2],
                                                    energy_mat[freq,3,subject_num-1]+metrics_save[3]]
                                    row = [subjects[subject_num-1],group_label[group],frequencylabels[freq],
                                            metrics_save[0]/2,metrics_save[1]/2,metrics_save[2]/2,metrics_save[3]/2,metrics_save[freq]/2]
                                    with open (file_metrics,'a') as csvfile:
                                        testwriter = csv.writer(csvfile,delimiter=',')
                                        testwriter.writerow(row)

                            row = [subjects[subject_num-1],group_label[group],frequencylabels[freq],
                                    energy_mat[freq,0,subject_num-1],energy_mat[freq,1,subject_num-1],
                                    energy_mat[freq,2,subject_num-1],energy_mat[freq,3,subject_num-1],
                                    energy_mat[freq,freq,subject_num-1]]

                            with open (file_metrics,'a') as csvfile:
                                testwriter = csv.writer(csvfile,delimiter=',')
                                testwriter.writerow(row)

        if make_plots==1:
            # Take average of each signal and re-normalizes
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

            if group == 0:
                freq = 0
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

            for freq in range(0,4):
                if group == 0:
                    mag_list.append(A_Fmag[0,:])
                    xy_xlist.append(A_Fx[0,:])
                    xy_ylist.append(A_Fy[0,:])

            ####################################################################
            ####################################################################
            ######   Plot freq spectrum for force-condition group   ############
            ####################################################################
            ####################################################################
            title = 'Frequency Spectrum Plot: ' + forceconditions[group]
            savename = 'Plots/'+'xy_freq_'+group_label[group]+'.png'
            if make_plot_each_sub==1:
                title = 'Sub' + str(sub_plot) + ' Frequency Spectrum Plot:\n' + forceconditions[group]
                savename = 'Plots/IndividualSubjectPlots/Sub'+str(sub_plot)+'_xy_freq_'+group_label[group]+'.png'
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
            fig.savefig(savename)

            title = 'Frequency Spectrum Plot: ' + forceconditions[group]
            savename = 'Plots/'+'mag_freq_'+group_label[group]+'.png'
            if make_plot_each_sub==1:
                title = 'Sub' + str(sub_plot) + ' Frequency Spectrum Plot:\n' + forceconditions[group]
                savename = 'Plots/IndividualSubjectPlots/Sub'+str(sub_plot)+'_mag_freq_'+group_label[group]+'.png'
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
            [fig,ax] = mag_spectrum(w,mag_list[n-4:n],freq_pendulum,title,xlabel,ylabel,legend,linestyles,colors,ymin,ymax,ymax_pend)
            fig.savefig(savename)
            ####################################################################
            ####################################################################
            ## Boxplot comparing frequency metrics for force-condition group ###
            ####################################################################
            ####################################################################
            if group != 0 and make_plot_each_sub==0:
                energy_mat = np.divide(energy_mat,energy_num)
                data = []
                for freq2 in range(0,4):
                    for freq in range(0,4):
                        energy = energy_mat[freq,freq2,:]
                        data.append(np.array(energy))

                # create lists for colors, labels, ect.
                # add lists for future plots
                n = len(data)
                labels = []
                box_alpha = []
                box_colors = []
                for i in range(n):
                    labels.append(freq_pendulum[i % 4])
                    box_alpha.append(1)
                    if (i // 4) == (i % 4): # if the ball and range match
                        box_colors.append('#ff7e0d')
                        # if group >= 2:
                        #     compare_groups1[group-2,(i//4),:] = data[i]
                        if group_label[group] == 'F1_B0':
                            compare_groups2[0,(i//4),:] = data[i]
                        elif group_label[group] == 'F0_B1':
                            compare_groups1[0,(i//4),:] = data[i]
                        elif group_label[group] == 'F1_B1':
                            compare_groups2[1,(i//4),:] = data[i]
                            compare_groups1[1,(i//4),:] = data[i]
                        compare_groups3[group,(i//4),:] = data[i]
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
        if make_plot_each_sub==0:
            for plot in range(3):
                figure_size = (5,3) # sets the size of the figure in inches
                xlabel = ''
                ylabel = 'Fraction of Total Energy'

                if plot==0:
                    title = 'Energy Content of Movement at Resonance\nFor Moving Ball Trials'
                    compare_groups = compare_groups1
                    legend_cond = ['Without Haptic Forces','With Haptic Forces']
                    # print(colors[1])
                    # colors = ['#601A4A','#EE442F','#63ACBE','#006400']
                    colors = ['#601A4A','#63ACBE']
                elif plot==1:
                    title = 'Energy Content of Movement at Resonance\nFor Haptic Force Trials'
                    compare_groups = compare_groups2
                    legend_cond = ['Still Ball','Moving Ball']
                    colors = ['#EE442F','#63ACBE']
                else:
                    title = 'Energy Content of Movement For All Experimental Conditions'
                    compare_groups = compare_groups3
                    legend_cond = group_label
                    colors = ['#601A4A','#EE442F','#63ACBE','#006400']
                n_group = compare_groups.shape[0]
                data = []
                labels = []
                box_alpha = []
                box_colors = []
                for freq in range(0,4):
                    for group in range(0,n_group):
                        data.append(np.array(compare_groups[group,freq,:]))
                        labels.append('')
                        box_alpha.append(1)
                        box_colors.append(colors[group])
                [fig,ax]=make_boxplot(data,title,xlabel,ylabel,labels,box_colors,box_alpha,figure_size)

                # Add rectangles for background
                [ymin,ymax]=ax.get_ylim()
                grey_colors = ['#e3e3e3','#ffffff','#e3e3e3','#ffffff']
                x1 = []
                x2 = []
                for i in range(4):
                    x = i*n_group+0.5
                    ax.add_patch(Rectangle((x,ymin), n_group, ymax-ymin, facecolor=grey_colors[i], alpha=1,zorder=1))
                    x1.append(x)
                    if i!=0:
                        x2.append(x)
                x2.append(x+n_group)

                # Text Labels
                y = ymin - ((ymax-ymin)/15)
                text_buffer = ((ymax-ymin)/12)
                name =frequencylabels
                add_labels(ax,x1,x2,y,name,text_buffer)
                y = ymin - ((ymax-ymin)/7.5)
                ax.text(0.5, y,'Ball:',horizontalalignment='right', fontname="Arial", fontsize=10, fontweight='bold')

                x = 5
                l1 = ax.add_patch(Rectangle((x,ymin-100), 3, ymax-ymin, facecolor=colors[0], alpha=1,zorder=1))
                l2 = ax.add_patch(Rectangle((x,ymin-100), 3, ymax-ymin, facecolor=colors[1], alpha=1,zorder=1))
                # l3 = ax.add_patch(Rectangle((x,ymin-100), 3, ymax-ymin, facecolor=colors[2], alpha=1,zorder=1))
                # l4 = ax.add_patch(Rectangle((x,ymin-100), 3, ymax-ymin, facecolor=colors[3], alpha=1,zorder=1))
                ax.set_ylim(ymin,ymax)

                if plot!=2:
                    # fig.subplots_adjust(right=0.72)
                    fig.subplots_adjust(bottom=0.28)
                    L = fig.legend([l1,l2], legend_cond,ncol=1, fontsize=10,loc='lower center')

                fig.savefig('Plots/'+'energy_comparegroups'+str(plot)+'.png')

        ####################################################################
        ####################################################################
        ######   Frequency spectrum plot for each frequency  #####
        ####################################################################
        ####################################################################
        title = 'Frequency Spectrum Plot'
        xlabel = 'Frequency (Hz)'
        ylabel = 'Normalized Frequency Amplitude'
        colors = ['#601A4A','#EE442F','#63ACBE','#006400']
        # colors = ['#EE442F','#63ACBE']
        linestyles = ['-','-','-','-']
        if haptic_forces_added:
            colors = ['#601A4A','#EE442F','#63ACBE','#006400','#ff7e0d']
            linestyles = ['-','-','-','-','-']
        ymin = 0.1
        ymax = 2
        ymax_pend = 1.5
        for freq in range(0,4):
            title = 'Frequency Spectrum Plot for ' + frequencylabels[freq]
            if make_plot_each_sub==1:
                title = 'Sub'+str(sub_plot)+' Frequency Spectrum Plot for ' + frequencylabels[freq]
            if haptic_forces_added:
                legend = ["Ball's resonant\nfrequency",
                            'Baseline Motion',
                            'Still Ball With\nHaptic Forces',
                            'Moving Ball Without\nHaptic Forces',
                            'Haptic Forces\nAdded to Signal',
                            'Moving Ball With\nHaptic Forces']
            else:
                legend = ["Ball's resonant\nfrequency",
                            'Baseline Motion',
                            'Still Ball With\nHaptic Forces',
                            'Moving Ball Without\nHaptic Forces',
                            'Moving Ball With\nHaptic Forces']
                # legend = ["Ball's resonant\nfrequency",
                #             'Still Ball With\nHaptic Forces',
                #             'Moving Ball With\nHaptic Forces']
            xdata = []
            ydata = []
            xdata.append(xy_xlist[freq])
            xdata.append(xy_xlist[freq+4])
            xdata.append(xy_xlist[freq+8])
            xdata.append(xy_xlist[freq+12])
            if haptic_forces_added:
                xdata.append(xy_xlist[freq+16])
            ydata.append(xy_ylist[freq])
            ydata.append(xy_ylist[freq+4])
            ydata.append(xy_ylist[freq+8])
            ydata.append(xy_ylist[freq+12])
            if haptic_forces_added:
                ydata.append(xy_ylist[freq+16])
            [fig,ax] = xy_spectrum(w,xdata,ydata,[freq_pendulum[freq]],title,xlabel,ylabel,legend,linestyles,colors,ymin,ymax,ymax_pend)
            savename = 'Plots/'+'xy_freq_'+frequencylabels[freq]+'.png'
            if make_plot_each_sub==1:
                savename = 'Plots/IndividualSubjectPlots/Sub'+str(sub_plot)+'_xy_freq_'+frequencylabels[freq]+'.png'
            fig.savefig(savename)

            data = []
            data.append(mag_list[freq])
            data.append(mag_list[freq+4])
            data.append(mag_list[freq+8])
            data.append(mag_list[freq+12])
            if haptic_forces_added:
                data.append(mag_list[freq+16])
            n = len(xy_xlist)
            [fig,ax] = mag_spectrum(w,data,[freq_pendulum[freq]],title,xlabel,ylabel,legend,linestyles,colors,ymin,ymax,ymax_pend)
            savename = 'Plots/'+'mag_freq_'+frequencylabels[freq]+'.png'
            if make_plot_each_sub==1:
                savename = 'Plots/IndividualSubjectPlots/Sub'+str(sub_plot)+'_mag_freq_'+frequencylabels[freq]+'.png'
            fig.savefig(savename)
# plt.show()
