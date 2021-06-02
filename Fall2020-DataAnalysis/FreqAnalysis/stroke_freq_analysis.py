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
# This program analyses the frequency spectrum for each trial/participant.
# This code expects the original data to be formatted and named properly
################################################################################

# Edit these variables before running
save_values = 1 # 0-do not save values for statistical tests 1-do save values
make_plots = 1 # 0-do not make plots 1-make plots
make_plot_each_sub = 0 # 0-do not make plots 1-make plots
filter = 1
window = .3
only_modsevere = 0


# subjects = [202,203,205,207,208,209,211,212]
subjects = [202,203,208,209,211,212,214]
subjects = [202,203,208,209,211,212]
make_line_plot = 0
if only_modsevere:
    subjects = [208,211,212]
    make_line_plot = 1
# subjects = [214,215]

DIR = "/home/mschlafly/Desktop/Stroke/" #set directory where data is mounted Ola- "/media/ola/Elements/R01prelim" Milli -"Z:"
DT = 0.05
Fs = 1/DT
freq_pendulum = [.5,1,1.5,2.5]
# s103flag = 0 # 0: 10% loading, 1: tabletop

# Label factors as strings for ezANOVA analysis
ind = [0,1,2,3]
frequencylabels = ['0.5Hz','1Hz','1.5Hz','2.5Hz']
markerstyles = ['o','D']
arms = ['paretic','nonparetic']
arms_label = ['A0','A1']
support_levels = ['0%','35%']
SL_label = ['SL1','SL2']
conditions = ['paretic-0%','paretic-35%','nonparetic-0%','nonparetic-35%']

colors = ['#998ec3', '#f1a340']
linestyles = ['--','-']


# Store data for statistical analysis
if save_values==1:
    file_metrics = "stroke-freq.csv"
    columns = ["Subject","Arm","Loading","BallFreq","EResonance"]
    with open(file_metrics,'w') as csvfile:
        testwriter = csv.writer(csvfile,delimiter=',')
        testwriter.writerow(columns)

    file_metrics_selected = "stroke-freq-average.csv"
    columns = ["Subject","BallFreq","SL1","SL0","A0","A1"]
    with open(file_metrics_selected,'w') as csvfile:
        testwriter = csv.writer(csvfile,delimiter=',')
        testwriter.writerow(columns)

    file_metrics_all = "stroke-freq-all.csv"
    columns = ["Subject","BallFreq","SL1","SL0","A0","A1"]
    with open(file_metrics_all,'w') as csvfile:
        testwriter = csv.writer(csvfile,delimiter=',')
        testwriter.writerow(columns)

# Create vector of frequencies of interest
nyquist_freq = int(np.floor(Fs/2))
print('Highest frequency evaluated: ', nyquist_freq)
nyquist_freq = 5
freq_step = 0.15 # set the resolution for the frequency bins

# Fill a vector of frequencies according to these specifications
w = np.zeros(int(np.floor(nyquist_freq/freq_step)-1))
frq = 0 # freq_step
for i in range(len(w)):
    w[i] = frq
    frq = frq + freq_step
print('w: ', w)
w_len = len(w)

def search_directory(DIR,sub_str,arm_str,SL_str):
    pathName = DIR + "S" + sub_str + "/"
    Files = []
    fileNames = os.listdir(pathName)
    for fileName in fileNames:
        if fileName.startswith("OutputData"):
            if "_Freq"+str(freq) in fileName:
                if arm_str in fileName:
                    if SL_str in fileName:
                            Files.append(fileName)
    return [pathName,Files]

# remove rows where the ball is green
def remove_lowe(df,freq):
    gravity = 9.81
    mass = 1.0
    percent_height = 0.3
    radius_options = [0.995, 0.249, 0.111, 0.04]
    R = radius_options[freq]
    max_energy = mass * gravity * percent_height * R
    # print('max_energy',max_energy)
    df = df[df['ball_energy'] > max_energy]
    return df

# remove rows where the person is not lifted
def remove_notlifted(df):
    table_z = -0.14
    z_tolerance = table_z + 0.01
    df = df[df['z'] > z_tolerance]
    return df


# if plotting each subject, iterate through each subject plot
num_sub_plots = 1
if make_plot_each_sub:
    num_sub_plots = len(subjects)
for sub_plot in range(0,num_sub_plots):

    # fill arrays for plotting the percentage decrease in function
    SL1 = np.zeros((len(frequencylabels),len(subjects)))
    A0 = np.zeros((len(frequencylabels),len(subjects)))
    SL0 = np.zeros((len(frequencylabels),len(subjects)))
    A1 = np.zeros((len(frequencylabels),len(subjects)))
    mean_energy = np.zeros((len(arms),len(support_levels),len(frequencylabels),len(subjects)))
    num_energy = np.zeros((len(arms),len(support_levels),len(frequencylabels),len(subjects)))

    if make_plots:

        # For aggregate frequency spectrum plot
        # A_Fmag_all = np.zeros((len(arms),len(support_levels),len(frequencylabels),len(subjects),w_len))
        A_Fmag = np.zeros((len(arms),len(support_levels),len(frequencylabels),len(subjects),w_len))
        # num_freq_all = np.zeros((len(arms),len(support_levels),len(frequencylabels),len(subjects)))
        num_freq = np.zeros((len(arms),len(support_levels),len(frequencylabels),len(subjects)))

    for subject_num in range(0,len(subjects)): #iterate though subjects

        if (make_plot_each_sub==0) or (make_plot_each_sub==1 and sub_plot==subject_num):

            # # For frequency spectrum plot
            # A_Fmag = np.zeros((len(arms),len(support_levels),len(frequencylabels),w_len))
            # num_freq = np.zeros((len(arms),len(support_levels),len(frequencylabels)))

            if make_plot_each_sub and make_plots:


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

            energy_array = np.zeros((len(arms),len(support_levels),len(frequencylabels),12))
            for arm in range(0,len(arms)):

                for SL in range(0,len(support_levels)):

                    arm_SL_list = []

                    for freq in range(0,len(frequencylabels)):

                        # if subjects[subject_num]==205 or subjects[subject_num]==207:
                        #     if freq==2 or freq==3:
                        #         print('skipping trial')
                        #         continue

                        energy_at_resonance = np.array([])

                        [pathName,Files] = search_directory(DIR,str(subjects[subject_num]),arms_label[arm],SL_label[SL])
                        trialnum_list = []
                        for i in Files:
                            str_split = i.split("_")
                            trial_num = int(str_split[2].split("Trial")[1])
                            trialnum_list.append(trial_num)
                        zipped_lists = zip(trialnum_list, Files)
                        sorted_zipped_lists = sorted(zipped_lists)

                        i_ee = 0
                        for i in range(len(sorted_zipped_lists)):

                            print(sorted_zipped_lists[i][1])
                            file = open(os.path.join(pathName, sorted_zipped_lists[i][1]), "rU")
                            data = genfromtxt(file,delimiter=',',dtype=float)
                            data = np.delete(data,0,0) # Deletes the first column of column names
                            df = pd.DataFrame({'z':data[:,3],'fx':data[:,8],'fy':data[:,9],'ball_energy':data[:,16]})
                            df = remove_lowe(df,freq)
                            df = remove_notlifted(df)
                            if len(df) > w_len:
                                y_Fx = df['fx'].tolist()
                                A_Fx_i,frq = calculate_amplitude(w,y_Fx,Fs)
                                y_Fy = df['fy'].tolist()
                                A_Fy_i,frq = calculate_amplitude(w,y_Fy,Fs)

                                # if freq == 2 and SL==0 and arm==0:
                                #     print(A_Fy_i)
                                # Normalize the spectrum so that the energy=1
                                dw = w[1]-w[0]
                                Ax_norm = normalize_spectrum(A_Fx_i,dw)
                                Ay_norm = normalize_spectrum(A_Fy_i,dw)

                                # Add the x and y spectrums together and normalize
                                Amag_norm = normalize_spectrum(Ax_norm+Ay_norm,dw)

                                # Check to see if the max amplitude is within 0.5 of desired amplitude
                                # Check to see if the energy@resonance > energy@other_frequencies
                                # use_trial = 1
                                # if freq!=0: # don't do this for 0.5Hz
                                #     window_keep = 0.5
                                #     max_index = np.argmax(Amag_norm)
                                #
                                #     # if w[max_index]>freq_pendulum[freq]+window_keep or w[max_index]<freq_pendulum[freq]-window_keep:
                                #     if abs(w[max_index]-freq_pendulum[freq])>window_keep:
                                #         # print('skip trial',w[max_index])
                                #         use_trial = 0

                                use_trial = 1
                                if use_trial:
                                    # Store the signal
                                    A_Fmag[arm,SL,freq,subject_num,:] += Amag_norm
                                    num_freq[arm,SL,freq,subject_num] += 1

                                    ee = find_energy_at_resonance(w,Amag_norm,freq_pendulum[freq],window)

                                    energy_at_resonance = np.append(energy_at_resonance,ee)
                                    energy_array[arm,SL,freq,i_ee] = ee
                                    i_ee += 1

                                    mean_energy[arm,SL,freq,subject_num] += ee
                                    num_energy[arm,SL,freq,subject_num] += 1

                                    if save_values:
                                        row = ["Sub"+str(subjects[subject_num]),arms[arm],support_levels[SL],frequencylabels[freq],ee]
                                        with open(file_metrics,'a') as csvfile:
                                            testwriter = csv.writer(csvfile,delimiter=',')
                                            testwriter.writerow(row)

                        arm_SL_list.append(energy_at_resonance)

                        # add subject's spectrum to aggregate plot
                        # A_Fmag_all[arm,SL,freq,:] += A_Fmag[arm,SL,freq,:]
                        # num_freq_all[arm,SL,freq] += 1

                    if make_plot_each_sub and make_plots:

                        # Plot energy@resonance scatterplot
                        data_mean = []
                        data_std = []
                        for freq in range(0,4):
                            if len(arm_SL_list)>3:
                                data_mean.append(np.mean(arm_SL_list[freq]))
                                data_std.append(np.std(arm_SL_list[freq])/np.sqrt(len(arm_SL_list[freq])))
                                # data_std.append(np.std(arm_SL_list[freq]))
                            else:
                                data_mean.append(0)
                                data_std.append(0)

                        i = arm *2 + SL
                        p[i] = ax_sub.errorbar(ind,data_mean,yerr=data_std,color=colors[arm],ls=linestyles[SL],marker="o",ecolor=colors[arm],capsize=5)

            # Normalize subject spectrums
            for freq in range(len(frequencylabels)):
                for arm in range(0,len(arms)):
                    for SL in range(0,len(support_levels)):
                        A_Fmag[arm,SL,freq,subject_num,:] /= num_freq[arm,SL,freq,subject_num]
                        A_Fmag[arm,SL,freq,subject_num,:] = normalize_spectrum(A_Fmag[arm,SL,freq,subject_num,:],dw)
                        # print("Energy: ",np.sum(np.square(A_Fmag[arm,SL,freq,subject_num,:])*dw))


            if save_values:
                # go through energy_array matching up values for percent decrease data
                for freq in range(0,len(frequencylabels)):
                    for i in range(12):
                        # np.zeros((len(arms),len(support_levels),len(frequencylabels),9))

                        # SL1
                        if min(abs(energy_array[:,1,freq,i]))>0.0:
                            paretic_ee = energy_array[0,1,freq,i]
                            nonparetic_ee = energy_array[1,1,freq,i]
                            SL1_diff = (nonparetic_ee - paretic_ee)/nonparetic_ee
                        else:
                            SL1_diff = 1000

                        # SL0
                        if min(abs(energy_array[:,0,freq,i]))>0.0:
                            paretic_ee = energy_array[0,0,freq,i]
                            nonparetic_ee = energy_array[1,0,freq,i]
                            SL0_diff = (nonparetic_ee - paretic_ee)/nonparetic_ee
                        else:
                            SL0_diff = 1000

                        # A0
                        if min(abs(energy_array[0,:,freq,i]))>0.0:
                            SL1_ee = energy_array[0,1,freq,i]
                            SL0_ee = energy_array[0,0,freq,i]
                            A0_diff = (SL0_ee - SL1_ee)/SL0_ee
                        else:
                            A0_diff = 1000

                        # A1
                        if min(abs(energy_array[1,:,freq,i]))>0.0:
                            SL1_ee = energy_array[1,1,freq,i]
                            SL0_ee = energy_array[1,0,freq,i]
                            A1_diff = (SL0_ee - SL1_ee)/SL0_ee
                        else:
                            A1_diff = 1000

                        row = ["Sub"+str(subjects[subject_num]),frequencylabels[freq],SL1_diff,SL0_diff,A0_diff,A1_diff]
                        with open(file_metrics_all,'a') as csvfile:
                            testwriter = csv.writer(csvfile,delimiter=',')
                            testwriter.writerow(row)

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
                ymin = 0
                ymax = 1
                ymax_pend = 1
                for freq in range(len(frequencylabels)):
                    mag_list = []
                    for arm in range(0,len(arms)):
                        for SL in range(0,len(support_levels)):
                            # A_Fmag[arm,SL,freq,subject_num,:] /= num_freq[arm,SL,freq,subject_num]
                            cutoff = 5  # desired cutoff frequency of the filter, Hz
                            # A_Fmag[arm,SL,freq,subject_num,:] = butter_lowpass_filter(A_Fmag[arm,SL,freq,subject_num,:], cutoff, Fs)
                            A_Fmag_temp= butter_lowpass_filter(A_Fmag[arm,SL,freq,subject_num,:], cutoff, Fs)
                            A_Fmag_temp= normalize_spectrum(A_Fmag_temp,dw)
                            # A_Fmag[arm,SL,freq,subject_num,:] = normalize_spectrum(A_Fmag[arm,SL,freq,subject_num,:],dw)
                            mag_list.append(A_Fmag_temp)
                    freq_plot = [freq_pendulum[freq]]
                    title = 'Force Frequency Spectrum for the '+frequencylabels[freq]+' Ball'
                    [fig_spec,ax_spec] = mag_spectrum(w,mag_list,freq_plot,title,xlabel,ylabel,legend,linestyles_spec,colors_spec,ymin,ymax,ymax_pend)
                    fig_spec.savefig('Plots/IndividualSubjectPlots/'+'S'+str(subjects[subject_num])+'/'+'S'+str(subjects[subject_num])+'_'+frequencylabels[freq]+'.pdf')


            if make_plots:
                for freq in range(0,len(frequencylabels)):
                    # print(num_freq[:,1,freq])
                    min_trials = 0
                    if min(num_freq[:,1,freq,subject_num])>=min_trials:
                        # store percentage loss in function for arm 35% loading
                        paretic_ee = find_energy_at_resonance(w,A_Fmag[0,1,freq,subject_num],freq_pendulum[freq],window)
                        nonparetic_ee = find_energy_at_resonance(w,A_Fmag[1,1,freq,subject_num],freq_pendulum[freq],window)

                        SL1[freq,subject_num] = (nonparetic_ee - paretic_ee)/nonparetic_ee
                    if min(num_freq[0,:,freq,subject_num])>=min_trials:
                        # store percentage loss in function with loading for paretic arm
                        SL1_ee = find_energy_at_resonance(w,A_Fmag[0,1,freq,subject_num],freq_pendulum[freq],window)
                        SL0_ee = find_energy_at_resonance(w,A_Fmag[0,0,freq,subject_num],freq_pendulum[freq],window)

                        A0[freq,subject_num] = (SL0_ee - SL1_ee)/SL0_ee
                    if min(num_freq[:,0,freq,subject_num])>=min_trials:
                        # store percentage loss in function for arm 0% loading
                        paretic_ee = find_energy_at_resonance(w,A_Fmag[0,0,freq,subject_num],freq_pendulum[freq],window)
                        nonparetic_ee = find_energy_at_resonance(w,A_Fmag[1,0,freq,subject_num],freq_pendulum[freq],window)

                        SL0[freq,subject_num] = (nonparetic_ee - paretic_ee)/nonparetic_ee

                    if min(num_freq[1,:,freq,subject_num])>=min_trials:
                        # store percentage loss in function with loading for nonparetic arm
                        SL1_ee = find_energy_at_resonance(w,A_Fmag[1,1,freq,subject_num],freq_pendulum[freq],window)
                        SL0_ee = find_energy_at_resonance(w,A_Fmag[1,0,freq,subject_num],freq_pendulum[freq],window)

                        A1[freq,subject_num] = (SL0_ee - SL1_ee)/SL0_ee

                    if save_values:
                        row = ["Sub"+str(subjects[subject_num]),frequencylabels[freq],SL1[freq,subject_num],SL0[freq,subject_num],A0[freq,subject_num],A1[freq,subject_num]]
                        with open(file_metrics_selected,'a') as csvfile:
                            testwriter = csv.writer(csvfile,delimiter=',')
                            testwriter.writerow(row)

    if make_plots and make_plot_each_sub==0:
        starting_i = 1


        # Aggregate frequency plot
        xlabel = ''
        ylabel = 'Frequency Amplitude'
        linestyles_spec = ['-','-','-','-']
        colors_spec = ['#f5793a','#a95aa1','#85c0f9','#0f2080']
        legend = ['Resonant Frequency','paretic-0%','paretic-35%','nonparetic-0%','nonparetic-35%']
        ymin = 0
        ymax = 1
        ymax_pend = 1
        for freq in range(len(frequencylabels)):
            mag_list = []
            for arm in range(0,len(arms)):
                for SL in range(0,len(support_levels)):
                    A_agg = np.zeros(w_len)
                    for subject_num in range(0,len(subjects)): #iterate though subjects
                        A_agg += A_Fmag[arm,SL,freq,subject_num,:]
                    A_agg /= len(subjects)
                    cutoff = 5  # desired cutoff frequency of the filter, Hz
                    A_agg = butter_lowpass_filter(A_agg, cutoff, Fs)
                    A_agg = normalize_spectrum(A_agg,dw)
                    mag_list.append(A_agg)

                    # A_Fmag_all[arm,SL,freq,:] /= num_freq_all[arm,SL,freq]
                    # cutoff = 5  # desired cutoff frequency of the filter, Hz
                    # A_Fmag_all[arm,SL,freq,:] = butter_lowpass_filter(A_Fmag_all[arm,SL,freq,:], cutoff, Fs)
                    # A_Fmag_all[arm,SL,freq,:] = normalize_spectrum(A_Fmag_all[arm,SL,freq,:],dw)
                    # mag_list.append(A_Fmag_all[arm,SL,freq,:])
            freq_plot = [freq_pendulum[freq]]
            title = 'Force Frequency Spectrum for the '+frequencylabels[freq]+' Ball'
            [fig_spec,ax_spec] = mag_spectrum(w,mag_list,freq_plot,title,xlabel,ylabel,legend,linestyles_spec,colors_spec,ymin,ymax,ymax_pend)
            if only_modsevere:
                fig_spec.savefig('Plots/agg_spectrum_'+frequencylabels[freq]+'_modsevere.pdf')

            else:
                fig_spec.savefig('Plots/agg_spectrum_'+frequencylabels[freq]+'.pdf')


    # Line plots for raw data
        figure_size = (6,3.55)
        fig_scat, ax_scat = plt.subplots(figsize=figure_size,dpi=300)
        xlabel = 'Resonant Frequency of Ball'
        ylabel = 'Fraction of Total Energy'
        title = 'Energy Content of Movement at Different Frequencies'

        # place grid in back
        ax_scat.grid(True, linestyle='-', which='major', axis='y', color='lightgrey',
                       alpha=0.5)
        ax_scat.set_axisbelow(True)

        # Add titles and labels
        plt.xlabel(xlabel,fontname="sans-serif", fontsize=11)
        plt.ylabel(ylabel,fontname="sans-serif", fontsize=11)
        plt.title(title,fontname="sans-serif", fontsize=11,fontweight='bold')
        for label in (ax_scat.get_yticklabels()):
            label.set_fontsize(8)

        # x-ticks x-axis
        plt.xticks(ind, frequencylabels, fontname="sans-serif", fontsize=10)
        for tick in ax_scat.get_xticklabels():
            tick.set_rotation(0)

        p = np.zeros(len(arms)+len(support_levels), dtype=object)

        mean_energy /= num_energy
        for arm in range(0,len(arms)):
            for SL in range(0,len(support_levels)):
                data_mean = []
                data_std = []
                for freq in range(0,4):
                    data_mean.append(np.mean(mean_energy[arm,SL,freq,:]))
                    data_std.append(np.std(mean_energy[arm,SL,freq,:])/np.sqrt(len(mean_energy[arm,SL,freq,:])))


                i = arm *2 + SL
                print(arm,SL,data_mean)
                p[i] = ax_scat.errorbar(ind,data_mean,yerr=data_std,color=colors[arm],ls=linestyles[SL],marker="o",ecolor=colors[arm],capsize=5)

        L = ax_scat.legend(p, conditions, fontsize=10,loc='upper right')
        # L = fig.legend(p, conditions, fontsize=10,loc='upper left',bbox_to_anchor=(0,0.35),bbox_transform=ax.transAxes)#'upper right')
        plt.setp(L.texts, family='sans-serif')
        fig_scat.subplots_adjust(bottom=0.15)
        if only_modsevere:
            fig_scat.savefig('Plots/e_at_res_raw_modsevere.pdf')
        else:
            fig_scat.savefig('Plots/e_at_res_raw.pdf')


        outlier = 3000000.
        # Make boxplot for SL1
        figure_size = (6,3.55) # sets the size of the figure in inches
        xlabel = 'Resonant Frequency of Ball'
        ylabel = 'Percentage Decrease in Energy (NP-P)/NP'
        title = 'Percentage Loss in Function in Paretic Arm With 35%Loading'
        box_colors = ['#FFFFFF','#FFFFFF','#FFFFFF','#FFFFFF']
        box_alpha = [1,1,1,1]
        labels = frequencylabels
        data = []
        for freq in range(0,len(frequencylabels)):
            freq_list = []
            for subject_num in range(0,len(subjects)):
                if abs(SL1[freq,subject_num])>0.00001 and abs(SL1[freq,subject_num])<outlier:
                    freq_list.append(SL1[freq,subject_num])
            data.append(np.array(freq_list))
        [fig,ax]=make_boxplot(data,title,xlabel,ylabel,labels,box_colors,box_alpha,figure_size)
        if make_plot_each_sub:
            fig.savefig('Plots/IndividualSubjectPlots/'+'S'+str(subjects[sub_plot])+'/'+'S'+str(subjects[sub_plot])+'_pl_SL1.pdf')
        else:
            fig.savefig('Plots/pl_SL1.pdf')
        if make_line_plot:
            # Line plot for SL1
            figure_size = (6,3.55) # sets the size of the figure in inches
            fig, ax = plt.subplots(figsize=figure_size,dpi=300)
            # place grid in back
            ax.grid(True, linestyle='-', which='major', axis='y', color='lightgrey',
                           alpha=0.5)
            ax.set_axisbelow(True)
            # Add titles and labels
            plt.xlabel(xlabel,fontname="sans-serif", fontsize=10)
            plt.ylabel(ylabel,fontname="sans-serif", fontsize=10)
            plt.title(title,fontname="sans-serif", fontsize=10,fontweight='bold')
            for label in (ax_scat.get_yticklabels()):
                label.set_fontsize(8)
            # x-ticks x-axis
            plt.xticks(ind, frequencylabels, fontname="sans-serif", fontsize=10)
            for tick in ax_scat.get_xticklabels():
                tick.set_rotation(0)

            p = np.zeros(len(subjects), dtype=object)

            for subject_num in range(0,len(subjects)):
                freq_list = []
                for freq in range(starting_i,len(frequencylabels)):
                    freq_list.append(SL1[freq,subject_num])
                p[subject_num] = ax.errorbar(ind[starting_i:4],freq_list)

            # if len(subjects)==7:
            #     sub_names = ['S202, FMA-51','S203, FMA-49','S208, FMA-37','S209, FMA-49','S211, FMA-17','S212, FMA-13','S214, FMA-30?']
            #     fig.subplots_adjust(right=0.75)
            #     L = fig.legend(p, sub_names, loc='center right', fontsize=9)
            #     plt.setp(L.texts, family='sans-serif')
            if only_modsevere:
                # sub_names = ['S208, FMA-37','S211, FMA-17','S212, FMA-13','S214, FMA-30?']
                sub_names = ['S208, FMA-37','S211, FMA-17','S212, FMA-13']
                fig.subplots_adjust(right=0.75)
                L = fig.legend(p, sub_names, loc='center right', fontsize=9)
                plt.setp(L.texts, family='sans-serif')
            fig.savefig('Plots/pl_SL1_modsevere.pdf')

        # Make boxplot for A0
        figure_size = (6,3.55) # sets the size of the figure in inches
        xlabel = 'Resonant Frequency of Ball'
        ylabel = 'Percentage Decrease in Energy (NL-L)/NL'
        title = 'Percentage Loss in Function From Loading in Paretic Arm'
        box_colors = ['#FFFFFF','#FFFFFF','#FFFFFF','#FFFFFF']
        box_alpha = [1,1,1,1]
        labels = frequencylabels
        data = []
        for freq in range(0,len(frequencylabels)):
            freq_list = []
            for subject_num in range(0,len(subjects)):
                if abs(A0[freq,subject_num])>0.00001 and abs(A0[freq,subject_num])<outlier:
                    freq_list.append(A0[freq,subject_num])
            data.append(np.array(freq_list))
        [fig,ax]=make_boxplot(data,title,xlabel,ylabel,labels,box_colors,box_alpha,figure_size)
        if make_plot_each_sub:
            fig.savefig('Plots/IndividualSubjectPlots/'+'S'+str(subjects[sub_plot])+'/'+'S'+str(subjects[sub_plot])+'_pl_A0.pdf')
        else:
            fig.savefig('Plots/pl_A0.pdf')


        if make_line_plot:
            # Line plot for A0
            figure_size = (6,3.55) # sets the size of the figure in inches
            fig, ax = plt.subplots(figsize=figure_size,dpi=300)
            # place grid in back
            ax.grid(True, linestyle='-', which='major', axis='y', color='lightgrey',
                           alpha=0.5)
            ax.set_axisbelow(True)
            # Add titles and labels
            plt.xlabel(xlabel,fontname="sans-serif", fontsize=10)
            plt.ylabel(ylabel,fontname="sans-serif", fontsize=10)
            plt.title(title,fontname="sans-serif", fontsize=10,fontweight='bold')
            for label in (ax_scat.get_yticklabels()):
                label.set_fontsize(8)
            # x-ticks x-axis
            plt.xticks(ind, frequencylabels, fontname="sans-serif", fontsize=10)
            for tick in ax_scat.get_xticklabels():
                tick.set_rotation(0)

            p = np.zeros(len(subjects), dtype=object)

            for subject_num in range(0,len(subjects)):
                freq_list = []
                for freq in range(starting_i,len(frequencylabels)):
                    freq_list.append(A0[freq,subject_num])
                p[subject_num] = ax.errorbar(ind[starting_i:4],freq_list)

            # if len(subjects)==7:
            #     sub_names = ['S202, FMA-51','S203, FMA-49','S208, FMA-37','S209, FMA-49','S211, FMA-17','S212, FMA-13','S214, FMA-30?']
            #     fig.subplots_adjust(right=0.75)
            #     L = fig.legend(p, sub_names, loc='center right', fontsize=9)
            #     plt.setp(L.texts, family='sans-serif')
            if only_modsevere:
                # sub_names = ['S208, FMA-37','S211, FMA-17','S212, FMA-13','S214, FMA-30?']
                sub_names = ['S208, FMA-37','S211, FMA-17','S212, FMA-13']
                fig.subplots_adjust(right=0.75)
                L = fig.legend(p, sub_names, loc='center right', fontsize=9)
                plt.setp(L.texts, family='sans-serif')
            fig.savefig('Plots/pl_A0_modsevere.pdf')



        # Make boxplot for SL0
        figure_size = (6,3.55) # sets the size of the figure in inches
        xlabel = 'Resonant Frequency of Ball'
        ylabel = 'Percentage Decrease in Energy (NP-P)/NP'
        title = 'Percentage Loss in Function in Paretic Arm With 0%Loading'
        box_colors = ['#FFFFFF','#FFFFFF','#FFFFFF','#FFFFFF']
        box_alpha = [1,1,1,1]
        labels = frequencylabels
        data = []
        for freq in range(0,len(frequencylabels)):
            freq_list = []
            for subject_num in range(0,len(subjects)):
                if abs(SL0[freq,subject_num])>0.00001 and abs(SL0[freq,subject_num])<outlier:
                    freq_list.append(SL0[freq,subject_num])
            data.append(np.array(freq_list))
        [fig,ax]=make_boxplot(data,title,xlabel,ylabel,labels,box_colors,box_alpha,figure_size)
        if make_plot_each_sub:
            fig.savefig('Plots/IndividualSubjectPlots/'+'S'+str(subjects[sub_plot])+'/'+'S'+str(subjects[sub_plot])+'_pl_SL0.pdf')
        else:
            fig.savefig('Plots/pl_SL0.pdf')



        if make_line_plot:
            # Line plot for SL0
            figure_size = (6,3.55) # sets the size of the figure in inches
            fig, ax = plt.subplots(figsize=figure_size,dpi=300)
            # place grid in back
            ax.grid(True, linestyle='-', which='major', axis='y', color='lightgrey',
                           alpha=0.5)
            ax.set_axisbelow(True)
            # Add titles and labels
            plt.xlabel(xlabel,fontname="sans-serif", fontsize=10)
            plt.ylabel(ylabel,fontname="sans-serif", fontsize=10)
            plt.title(title,fontname="sans-serif", fontsize=10,fontweight='bold')
            for label in (ax_scat.get_yticklabels()):
                label.set_fontsize(8)
            # x-ticks x-axis
            plt.xticks(ind, frequencylabels, fontname="sans-serif", fontsize=10)
            for tick in ax_scat.get_xticklabels():
                tick.set_rotation(0)

            p = np.zeros(len(subjects), dtype=object)

            for subject_num in range(0,len(subjects)):
                freq_list = []
                for freq in range(starting_i,len(frequencylabels)):
                    freq_list.append(SL0[freq,subject_num])
                p[subject_num] = ax.errorbar(ind[starting_i:4],freq_list)

            # if len(subjects)==7:
            #     sub_names = ['S202, FMA-51','S203, FMA-49','S208, FMA-37','S209, FMA-49','S211, FMA-17','S212, FMA-13','S214, FMA-30?']
            #     fig.subplots_adjust(right=0.75)
            #     L = fig.legend(p, sub_names, loc='center right', fontsize=9)
            #     plt.setp(L.texts, family='sans-serif')
            if only_modsevere:
                # sub_names = ['S208, FMA-37','S211, FMA-17','S212, FMA-13','S214, FMA-30?']
                sub_names = ['S208, FMA-37','S211, FMA-17','S212, FMA-13']
                fig.subplots_adjust(right=0.75)
                L = fig.legend(p, sub_names, loc='center right', fontsize=9)
                plt.setp(L.texts, family='sans-serif')
            fig.savefig('Plots/pl_SL0_modsevere.pdf')


        # Make boxplot for A1
        figure_size = (6,3.55) # sets the size of the figure in inches
        xlabel = 'Resonant Frequency of Ball'
        ylabel = 'Percentage Decrease in Energy (NL-L)/NL'
        title = 'Percentage Loss in Function From Loading in Non-Paretic Arm'
        box_colors = ['#FFFFFF','#FFFFFF','#FFFFFF','#FFFFFF']
        box_alpha = [1,1,1,1]
        labels = frequencylabels
        data = []
        for freq in range(0,len(frequencylabels)):
            freq_list = []
            for subject_num in range(0,len(subjects)):
                if abs(A1[freq,subject_num])>0.00001 and abs(A1[freq,subject_num])<outlier:
                    freq_list.append(A1[freq,subject_num])
            data.append(np.array(freq_list))
        [fig,ax]=make_boxplot(data,title,xlabel,ylabel,labels,box_colors,box_alpha,figure_size)
        if make_plot_each_sub:
            fig.savefig('Plots/IndividualSubjectPlots/'+'S'+str(subjects[sub_plot])+'/'+'S'+str(subjects[sub_plot])+'_pl_A1.pdf')
        else:
            fig.savefig('Plots/pl_A1.pdf')

        if make_line_plot:
            # Line plot for A1
            figure_size = (6,3.55) # sets the size of the figure in inches
            fig, ax = plt.subplots(figsize=figure_size,dpi=300)
            # place grid in back
            ax.grid(True, linestyle='-', which='major', axis='y', color='lightgrey',
                           alpha=0.5)
            ax.set_axisbelow(True)
            # Add titles and labels
            plt.xlabel(xlabel,fontname="sans-serif", fontsize=10)
            plt.ylabel(ylabel,fontname="sans-serif", fontsize=10)
            plt.title(title,fontname="sans-serif", fontsize=10,fontweight='bold')
            for label in (ax_scat.get_yticklabels()):
                label.set_fontsize(8)
            # x-ticks x-axis
            plt.xticks(ind, frequencylabels, fontname="sans-serif", fontsize=10)
            for tick in ax_scat.get_xticklabels():
                tick.set_rotation(0)

            p = np.zeros(len(subjects), dtype=object)

            for subject_num in range(0,len(subjects)):
                freq_list = []
                for freq in range(starting_i,len(frequencylabels)):
                    freq_list.append(A1[freq,subject_num])
                p[subject_num] = ax.errorbar(ind[starting_i:4],freq_list)

            # if len(subjects)==7:
            #     sub_names = ['S202, FMA-51','S203, FMA-49','S208, FMA-37','S209, FMA-49','S211, FMA-17','S212, FMA-13','S214, FMA-30?']
            #     fig.subplots_adjust(right=0.75)
            #     L = fig.legend(p, sub_names, loc='center right', fontsize=9)
            #     plt.setp(L.texts, family='sans-serif')
            if only_modsevere:
                # sub_names = ['S208, FMA-37','S211, FMA-17','S212, FMA-13','S214, FMA-30?']
                sub_names = ['S208, FMA-37','S211, FMA-17','S212, FMA-13']
                fig.subplots_adjust(right=0.75)
                L = fig.legend(p, sub_names, loc='center right', fontsize=9)
                plt.setp(L.texts, family='sans-serif')
            fig.savefig('Plots/pl_A1_modsevere.pdf')


print(SL1)
print(A0)
print(SL0)
print(A1)


# plt.show()
