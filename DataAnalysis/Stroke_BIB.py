import numpy as np
import scipy as sp
from numpy import fft, genfromtxt
from scipy import fft, arange, stats
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from utils.plot_utils import *
from utils.plot_freq_utils import *
from utils.perform_transform import *
import csv
import pandas as pd
import os

################################################################################
# This program analyses the frequency spectrum for each trial/participant.
# This code expects the original data to be formatted and named properly
################################################################################

# Edit these variables before running
save_values = 1 # 0-do not save values for statistical tests 1-do save values
make_plots = 1 # 0-do not make plots 1-make plots
make_plot_each_sub = 0 # 0-do not make plots 1-make plots; consider using subject list that contains excluded participants
if make_plot_each_sub:
    print('Note: Some aggregate plots will not be generated if make_plot_each_sub is True')
    if not make_plots:
        print('Individual plots will not generate unless make_plots is True')
window = .3000000000001
group = 'all' # 'all', mild stroke: 'm', moderate-severe: 'ms'
make_line_plot = 1
num_included_trials = 9 # number of trials included in statistics and plot per participant

# Included participants
subjects = [202,203,208,209,211,212,214,218,219,220]
FMA_list = [51,49,37,49,17,13,32,7,9,21]

# # With excluded participants
# subjects = [202,203,205,207,208,209,211,212,214,215,216,217,218,219,220]
# FMA_list = [51,49,34,30,37,49,17,13,32,0,0,0,7,9,21]

folder_name = 'all'
if group=='ms':
    folder_name = 'moderate-severe'
    new_sub_list = []
    new_FMA_list = []
    for i in range(len(subjects)):
        # if FMA_list[i]<35:
        if FMA_list[i]<40:
            new_sub_list.append(subjects[i])
            new_FMA_list.append(FMA_list[i])
    subjects = new_sub_list
    FMA_list = new_FMA_list
    save_values = 0
elif group=='m':
    folder_name = 'mild'
    new_sub_list = []
    new_FMA_list = []
    for i in range(len(subjects)):
        if FMA_list[i]>=40:
            new_sub_list.append(subjects[i])
            new_FMA_list.append(FMA_list[i])
    subjects = new_sub_list
    FMA_list = new_FMA_list
    save_values = 0

DIR = "/home/milli/Desktop/Stroke/" #set directory where data is mounted Ola- "/media/ola/Elements/R01prelim" Milli -"Z:"
DT = 0.05
Fs = 1/DT
freq_pendulum = [.5,1,1.5,2.5]
cutoff = 7  # desired cutoff frequency of the filter, Hz

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
if make_line_plot:
    sub_names = []
    for i in range(len(subjects)):
        # sub_names.append(str(subjects[i])+', FMA-'+str(FMA_list[i]))
        sub_names.append('FMA: '+str(FMA_list[i]))

# Create vector of frequencies of interest
nyquist_freq = int(np.floor(Fs/2))
print('Highest frequency evaluated: ', nyquist_freq)
nyquist_freq = 4
freq_step = 0.15 # set the resolution for the frequency bins

# Fill a vector of frequencies according to these specifications
w = np.zeros(int(np.floor(nyquist_freq/freq_step)-1))
frq = 0 # freq_step
for i in range(len(w)):
    w[i] = frq
    frq = frq + freq_step
print('w: ', w)
w_len = len(w)
dw = w[1]-w[0]

# Store data for statistical analysis
if save_values==1:
    file_metrics = "stroke-freq.csv"
    columns = ["Subject","FMA","Severity","Arm","Loading","BallFreq","EResonance"]
    with open(file_metrics,'w') as csvfile:
        testwriter = csv.writer(csvfile,delimiter=',')
        testwriter.writerow(columns)

    file_spectrum = "stroke-freq-spectrum.csv"
    columns = ["Subject","FMA","Arm","Loading","BallFreq"]
    for i in range(len(w)):
        columns.append(w[i])
    with open(file_spectrum,'w') as csvfile:
        testwriter = csv.writer(csvfile,delimiter=',')
        testwriter.writerow(columns)

    file_metrics_selected = "stroke-freq-percent-loss-aggregate.csv"
    columns = ["Subject","FMA","BallFreq","SL1","SL0","A0","A1","A_com"]
    with open(file_metrics_selected,'w') as csvfile:
        testwriter = csv.writer(csvfile,delimiter=',')
        testwriter.writerow(columns)

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

# if plotting each subject, iterate through each subject plot
num_sub_plots = 1
if make_plot_each_sub:
    num_sub_plots = len(subjects)
for sub_plot in range(0,num_sub_plots):

    # fill arrays for plotting the percentage decrease in function
    SL1 = np.zeros((len(frequencylabels),len(subjects))) # compare arms at SL1
    A0 = np.zeros((len(frequencylabels),len(subjects))) # compare SLs at A0
    SL0 = np.zeros((len(frequencylabels),len(subjects))) # compare arms at SL0
    A1 = np.zeros((len(frequencylabels),len(subjects))) # compare SLs at A1
    A_com = np.zeros((len(frequencylabels),len(subjects))) # compare arms with support levels combined
    mean_energy = np.zeros((len(arms),len(support_levels),len(frequencylabels),len(subjects)))
    num_energy = np.zeros((len(arms),len(support_levels),len(frequencylabels),len(subjects)))

    if make_plots:

        # For aggregate frequency spectrum plot
        A_Fmag = np.zeros((len(arms),len(support_levels),len(frequencylabels),len(subjects),w_len))
        num_freq = np.zeros((len(arms),len(support_levels),len(frequencylabels),len(subjects)))

    for subject_num in range(0,len(subjects)): #iterate though subjects
        if FMA_list[subject_num]<40:
            severity = 'mod-severe'
        else:
            severity = 'mild'

        if (make_plot_each_sub==0) or (make_plot_each_sub==1 and sub_plot==subject_num):

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
                plt.xlabel(xlabel,fontname="sans-serif", fontsize=9)
                plt.ylabel(ylabel,fontname="sans-serif", fontsize=9)
                plt.title(title,fontname="sans-serif", fontsize=9,fontweight='bold')
                for label in (ax_sub.get_yticklabels()):
                    label.set_fontsize(8)

                # x-ticks x-axis
                plt.xticks(ind, frequencylabels, fontname="sans-serif", fontsize=9)
                for tick in ax_sub.get_xticklabels():
                    tick.set_rotation(0)

                p = np.zeros(len(arms)+len(support_levels), dtype=object)

            energy_array = np.zeros((len(arms),len(support_levels),len(frequencylabels),num_included_trials))

            arm_list = [1,0]
            for arm in arm_list:

                for SL in range(0,len(support_levels)):

                    arm_SL_list = []

                    for freq in range(0,len(frequencylabels)):

                        energy_at_resonance = np.array([])

                        # Get all trials
                        [pathName,Files] = search_directory(DIR,str(subjects[subject_num]),arms_label[arm],SL_label[SL])
                        trialnum_list = []
                        for i in Files:
                            str_split = i.split("_")
                            trial_num = int(str_split[2].split("Trial")[1])
                            trialnum_list.append(trial_num)
                        zipped_lists = zip(trialnum_list, Files)
                        sorted_zipped_lists = sorted(zipped_lists)

                        # Decide which trials to include
                        sorted_zipped_lists_include = []
                        for i in range(len(sorted_zipped_lists)):
                            file = open(os.path.join(pathName, sorted_zipped_lists[i][1]), "r")
                            data = genfromtxt(file,delimiter=',',dtype=float)
                            data = np.delete(data,0,0) # Deletes the first column of column names
                            df = pd.DataFrame({'z':data[:,3],'fx':data[:,8],'fy':data[:,9],'ball_energy':data[:,16]})
                            df = remove_lowe(df,freq)
                            df = remove_notlifted(df)
                            if len(df) > w_len:
                                sorted_zipped_lists_include.append(sorted_zipped_lists[i])
                            else:
                                print('w\'s resolution is too high for ', sorted_zipped_lists[i][1])
                        final_included_trials = []
                        # print(len(sorted_zipped_lists_include))
                        if len(sorted_zipped_lists_include)<num_included_trials:
                            print('\n\n\n ERROR: not enough trials for this condition \n\n\n')
                        if len(sorted_zipped_lists_include)<=num_included_trials:
                            final_included_trials = sorted_zipped_lists_include
                        else:
                            for i in range(len(sorted_zipped_lists_include)):
                                # Skip trial/2, unless trial is missing
                                if len(sorted_zipped_lists_include)>=9: # full dataset
                                    if i<9-num_included_trials:
                                        continue
                                else:
                                    if i<len(sorted_zipped_lists_include)-num_included_trials:
                                        continue
                                final_included_trials.append(sorted_zipped_lists_include[i])
                                if len(final_included_trials)==num_included_trials:
                                    break

                        # print(final_included_trials)
                        i_ee = 0
                        for i in range(len(final_included_trials)):

                            print(final_included_trials[i][1])
                            file = open(os.path.join(pathName, final_included_trials[i][1]), "r")
                            data = genfromtxt(file,delimiter=',',dtype=float)
                            data = np.delete(data,0,0) # Deletes the first column of column names
                            df = pd.DataFrame({'z':data[:,3],'fx':data[:,8],'fy':data[:,9],'ball_energy':data[:,16]})
                            df = remove_lowe(df,freq)
                            df = remove_notlifted(df)

                            y_Fx = df['fx'].tolist()
                            A_Fx_i,frq = calculate_amplitude(w,y_Fx,Fs)
                            y_Fy = df['fy'].tolist()
                            A_Fy_i,frq = calculate_amplitude(w,y_Fy,Fs)

                            # Normalize the spectrum so that the energy=1
                            dw = w[1]-w[0]
                            Ax_norm = normalize_spectrum(A_Fx_i,dw)
                            Ay_norm = normalize_spectrum(A_Fy_i,dw)

                            # Add the x and y spectrums together and normalize
                            Amag_norm = normalize_spectrum(Ax_norm+Ay_norm,dw)

                            # if not window_plots:

                            if make_plots:
                                A_Fmag[arm,SL,freq,subject_num,:] += Amag_norm
                                num_freq[arm,SL,freq,subject_num] += 1

                            ee = find_energy_at_resonance(w,Amag_norm,freq_pendulum[freq],window)

                            energy_at_resonance = np.append(energy_at_resonance,ee)
                            energy_array[arm,SL,freq,i_ee] = ee
                            i_ee += 1

                            mean_energy[arm,SL,freq,subject_num] += ee
                            num_energy[arm,SL,freq,subject_num] += 1

                            if save_values:
                                row = ["Sub"+str(subjects[subject_num]),FMA_list[subject_num],severity,arms[arm],support_levels[SL],frequencylabels[freq],ee]
                                with open(file_metrics,'a') as csvfile:
                                    testwriter = csv.writer(csvfile,delimiter=',')
                                    testwriter.writerow(row)

                        arm_SL_list.append(energy_at_resonance)

                    if make_plot_each_sub and make_plots:

                        # Plot energy@resonance scatterplot
                        data_mean = []
                        data_std = []
                        for freq in range(0,4):
                            if len(arm_SL_list)>3:
                                data_mean.append(np.mean(arm_SL_list[freq]))
                                data_std.append(np.std(arm_SL_list[freq])/np.sqrt(len(arm_SL_list[freq])))
                            else:
                                print('Subject ', subjects[sub_plot], 'has less than 3 good trials')
                                data_mean.append(0)
                                data_std.append(0)

                        i = arm *2 + SL
                        p[i] = ax_sub.errorbar(ind,data_mean,yerr=data_std,color=colors[arm],ls=linestyles[SL],marker="o",ecolor=colors[arm],capsize=5)

            # Normalize subject spectrums
            # if not window_plots:
            for freq in range(len(frequencylabels)):
                for arm in range(0,len(arms)):
                    for SL in range(0,len(support_levels)):
                        A_Fmag[arm,SL,freq,subject_num,:] /= num_freq[arm,SL,freq,subject_num] # averaging first does not change the numbers at all
                        A_Fmag[arm,SL,freq,subject_num,:] = normalize_spectrum(A_Fmag[arm,SL,freq,subject_num,:],dw)
                        # print("Energy: ",np.sum(np.square(A_Fmag[arm,SL,freq,subject_num,:])*dw))

                        if save_values:
                            row = ["Sub"+str(subjects[subject_num]),FMA_list[subject_num],arms[arm],support_levels[SL],frequencylabels[freq]]
                            for i in range(len(w)):
                                row.append(A_Fmag[arm,SL,freq,subject_num,i])
                            with open(file_spectrum,'a') as csvfile:
                                testwriter = csv.writer(csvfile,delimiter=',')
                                testwriter.writerow(row)

            if make_plot_each_sub and make_plots:

                # Energy@resonance scatter plot add legend and save
                L = ax_sub.legend(p, conditions, fontsize=9,loc='upper right')
                # L = fig.legend(p, conditions, fontsize=9,loc='upper left',bbox_to_anchor=(0,0.35),bbox_transform=ax.transAxes)#'upper right')
                plt.setp(L.texts, family='sans-serif')
                # fig.subplots_adjust(right=0.7)
                fig_sub.savefig('Plots_stroke/IndividualSubjectPlots/'+'S'+str(subjects[subject_num])+'/'+'S'+str(subjects[subject_num])+'.pdf')
                plt.close(fig_sub)

                # Make frequency spectrum plots for every ball frequency
                xlabel = ''
                ylabel = 'Frequency Amplitude'
                linestyles_spec = ['-','-','-','-']
                colors_spec = ['#85c0f9','#0f2080','#f5793a','#a95aa1']
                legend = ['Resonant Frequency','nonparetic-0%','nonparetic-35%','paretic-0%','paretic-35%']
                ymin = 0
                ymax = 1.4
                ymax_pend = 1.4
                for freq in range(len(frequencylabels)):
                    mag_list = []
                    for arm in [1,0]:
                        for SL in range(0,len(support_levels)):
                            A_Fmag_temp= butter_lowpass_filter(A_Fmag[arm,SL,freq,subject_num,:], cutoff, Fs)
                            A_Fmag_temp= normalize_spectrum(A_Fmag_temp,dw)
                            mag_list.append(A_Fmag_temp)
                    freq_plot = [freq_pendulum[freq]]
                    title = 'Force Frequency Spectrum for the '+frequencylabels[freq]+' Ball'
                    [fig_spec,ax_spec] = mag_spectrum(w,mag_list,freq_plot,title,xlabel,ylabel,legend,linestyles_spec,colors_spec,ymin,ymax,ymax_pend)
                    fig_spec.savefig('Plots_stroke/IndividualSubjectPlots/'+'S'+str(subjects[subject_num])+'/'+'S'+str(subjects[subject_num])+'_'+frequencylabels[freq]+'.pdf')
                    ax_spec.set_xscale('log')
                    fig_spec.savefig('Plots_stroke/IndividualSubjectPlots/'+'S'+str(subjects[subject_num])+'/'+'S'+str(subjects[subject_num])+'_'+frequencylabels[freq]+'_log.pdf')
                    plt.close(fig_spec)

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

                    if min(num_freq[:,1,freq,subject_num])>=min_trials and min(num_freq[:,0,freq,subject_num])>=min_trials:
                        paretic_combined = normalize_spectrum(A_Fmag[0,0,freq,subject_num]+A_Fmag[0,1,freq,subject_num],dw)
                        paretic_ee = find_energy_at_resonance(w,paretic_combined,freq_pendulum[freq],window)

                        nonparetic_combined = normalize_spectrum(A_Fmag[1,0,freq,subject_num]+A_Fmag[1,1,freq,subject_num],dw)
                        nonparetic_ee = find_energy_at_resonance(w,nonparetic_combined,freq_pendulum[freq],window)

                        A_com[freq,subject_num] = (nonparetic_ee - paretic_ee)/nonparetic_ee

                    if save_values:
                        row = ["Sub"+str(subjects[subject_num]),FMA_list[subject_num],frequencylabels[freq],SL1[freq,subject_num],SL0[freq,subject_num],A0[freq,subject_num],A1[freq,subject_num],A_com[freq,subject_num]]
                        with open(file_metrics_selected,'a') as csvfile:
                            testwriter = csv.writer(csvfile,delimiter=',')
                            testwriter.writerow(row)

    if make_plots and make_plot_each_sub==0:
        starting_i = 1

        # Aggregate frequency plot
        xlabel = ''
        ylabel = 'Frequency Amplitude'
        linestyles_spec = ['-','-','-','-']
        colors_spec = ['#85c0f9','#0f2080','#f5793a','#a95aa1']
        legend = ['Resonant Frequency','nonparetic-0%','nonparetic-35%','paretic-0%','paretic-35%']
        ymin = 0
        ymax = 1.25
        ymax_pend = 1.25
        for freq in range(len(frequencylabels)):
            mag_list = []
            for arm in [1,0]:
                for SL in range(0,len(support_levels)):
                    A_agg = np.zeros(w_len)
                    for subject_num in range(0,len(subjects)): #iterate though subjects
                        A_agg += A_Fmag[arm,SL,freq,subject_num,:]
                        print(subject_num,A_Fmag[arm,SL,freq,subject_num,:])
                    A_agg /= len(subjects)
                    A_agg = butter_lowpass_filter(A_agg, cutoff, Fs)
                    A_agg = normalize_spectrum(A_agg,dw)
                    mag_list.append(A_agg)
            freq_plot = [freq_pendulum[freq]]
            title = 'Force Frequency Spectrum for the '+frequencylabels[freq]+' Ball'
            [fig_spec,ax_spec] = mag_spectrum(w,mag_list,freq_plot,title,xlabel,ylabel,legend,linestyles_spec,colors_spec,ymin,ymax,ymax_pend)
            fig_spec.savefig('Plots_stroke/'+folder_name+'/agg_spectrum_'+frequencylabels[freq]+'_'+group+'.pdf')
            fig_spec.savefig('Plots_stroke/'+folder_name+'/agg_spectrum_'+frequencylabels[freq]+'_'+group+'.png')
            ax_spec.set_xscale('log')
            fig_spec.savefig('Plots_stroke/'+folder_name+'/agg_spectrum_log_'+frequencylabels[freq]+'_'+group+'.pdf')
            fig_spec.savefig('Plots_stroke/'+folder_name+'/agg_spectrum_log_'+frequencylabels[freq]+'_'+group+'.png')
            plt.close(fig_spec)


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
        plt.xlabel(xlabel,fontname="sans-serif", fontsize=9)
        plt.ylabel(ylabel,fontname="sans-serif", fontsize=9)
        plt.title(title,fontname="sans-serif", fontsize=9,fontweight='bold')
        for label in (ax_scat.get_yticklabels()):
            label.set_fontsize(8)
        # x-ticks x-axis
        plt.xticks(ind, frequencylabels, fontname="sans-serif", fontsize=9)
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

        L = ax_scat.legend(p, conditions, fontsize=9,loc='upper right')
        # L = fig.legend(p, conditions, fontsize=9,loc='upper left',bbox_to_anchor=(0,0.35),bbox_transform=ax.transAxes)#'upper right')
        plt.setp(L.texts, family='sans-serif')
        fig_scat.subplots_adjust(bottom=0.15)
        fig_scat.savefig('Plots_stroke/'+folder_name+'/e_at_res_raw_'+group+'.pdf')
        fig_scat.savefig('Plots_stroke/'+folder_name+'/e_at_res_raw_'+group+'.png')
        plt.close(fig_scat)



        # Make boxplot for A_com
        figure_size = (3,3) # sets the size of the figure in inches
        xlabel = 'Resonant Frequency of Ball'
        ylabel = 'Percentage Decrease in Energy'
        title = 'Percentage Loss in Function\nin Paretic Arm'
        box_colors = ['#FFFFFF','#FFFFFF','#FFFFFF']
        box_alpha = [1,1,1]
        labels = frequencylabels[1:]
        data = []
        for freq in range(1,len(frequencylabels)): # skip 0.5Hz
            freq_list = []
            for subject_num in range(0,len(subjects)):
                # if abs(A_com[freq,subject_num])>0.00001 and abs(A_com[freq,subject_num])<outlier:
                freq_list.append(A_com[freq,subject_num])
            data.append(np.array(freq_list))
        [fig,ax] = make_boxplot(data,title,xlabel,ylabel,labels,box_colors,box_alpha,figure_size)
        fig.subplots_adjust(left=0.25)
        if group=='ms' and make_plot_each_sub==0:
            # Add stastical signicant marking #########################################
            sig_matrix = np.array([[1,2,.003],[0,2,0.017]])
            add_stats(data,sig_matrix,ax,spread_factor=20)
        fig.savefig('Plots_stroke/'+folder_name+'/pl_A_com_'+group+'.pdf')
        fig.savefig('Plots_stroke/'+folder_name+'/pl_A_com_'+group+'.png')
        plt.close(fig)
        if make_line_plot:
            figure_size = (5,3) # sets the size of the figure in inches
            fig, ax = plt.subplots(figsize=figure_size,dpi=300)
            # place grid in back
            ax.grid(True, linestyle='-', which='major', axis='y', color='lightgrey',
                           alpha=0.5)
            ax.set_axisbelow(True)
            # Add titles and labels
            plt.xlabel(xlabel,fontname="sans-serif", fontsize=9)
            plt.ylabel(ylabel,fontname="sans-serif", fontsize=9)
            plt.title(title,fontname="sans-serif", fontsize=9,fontweight='bold')
            for label in (ax_scat.get_yticklabels()):
                label.set_fontsize(8)
            # x-ticks x-axis
            plt.xticks(ind, frequencylabels, fontname="sans-serif", fontsize=9)
            for tick in ax_scat.get_xticklabels():
                tick.set_fontsize(8)
                tick.set_rotation(0)
            p = np.zeros(len(subjects), dtype=object)
            for subject_num in range(0,len(subjects)):
                freq_list = []
                for freq in range(starting_i,len(frequencylabels)):
                    freq_list.append(A_com[freq,subject_num])
                p[subject_num] = ax.errorbar(ind[starting_i:4],freq_list)
            fig.subplots_adjust(right=0.75)
            L = fig.legend(p, sub_names, loc='center right', fontsize=9)
            plt.setp(L.texts, family='sans-serif')
            fig.subplots_adjust(left=0.2,right=.75,bottom=0.2)
            fig.savefig('Plots_stroke/'+folder_name+'/pl_A_com_indiv_'+group+'.pdf')
            fig.savefig('Plots_stroke/'+folder_name+'/pl_A_com_indiv_'+group+'.png')
            plt.close(fig)

# plt.show()
