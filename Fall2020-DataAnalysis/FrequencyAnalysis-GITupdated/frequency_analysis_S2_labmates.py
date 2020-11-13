import numpy as np
import scipy as sp
from numpy import fft, genfromtxt
from scipy import fft, arange, stats
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from make_boxplot import make_boxplot, add_stats
from plot_freq_spectrum import xy_spectrum,mag_spectrum
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
make_plot_each_sub = 0 # 0-do not make plots 1-make plots
number_of_subjects = 7
DIR = "Z:Fall2020Data/" #set directory where data is mounted Ola- "/media/ola/Elements/R01prelim" Milli -"Z:"
DT = 0.05
Fs = 1/DT
freq_pendulum = [.5,1,1.5,2]
only_noforces = 1
split_signal = 1
freq_range_boxplot = 0

# Label factors as strings for ezANOVA analysis
hapticforces = ['off','on']
frequencylabels = ['0.5Hz','1Hz','1.5Hz','2Hz']
subjects = ['Subject1','Subject2','Subject3','Subject4','Subject5','Subject6','Subject7','Subject8','Subject9','Subject10','Subject11','Subject12','Subject13','Subject14','Subject15','Subject16','Subject17','Subject18']
# colors = ['#601A4A','#EE442F','#63ACBE','#006400'] # for plots
# linestyles = ['--','-']
markerstyles = ['o','D']
radius_options = [ 0.995, 0.249, 0.111, 0.062, 0.04 ];

# Create vector of frequencies of interest
nyquist_freq = int(np.floor(Fs/2))
# nyquist_freq = 40 # override the nyquist frequency (there is little info of value at higher frequencies)
print('Highest frequency evaluated: ', nyquist_freq)
freq_step = 0.15 # set the resolution for the frequcny bins

# Fill a vector of frequencies according to these specifications
w = np.zeros(int(np.floor(nyquist_freq/freq_step)-1))
frq = 0 # freq_step
for i in range(len(w)):
    w[i] = frq
    frq = frq + freq_step
print('w: ', w)
w_len = len(w)

# create matrices for saving energy
if freq_range_boxplot:
    window = .2 # for calculating the energy around different frequencies
    energy_mat_fx = np.zeros((4,4,6))
    energy_mat_fy = np.zeros((4,4,6))
    energy_num = np.zeros((4,4,6)) # saves the number of times trials for each condition

if make_plots:
    xy_xlist = []
    xy_ylist = []
    mag_list = []
# if plotting each subject, iterate through each subject plot
num_iters = 2
if make_plot_each_sub or freq_range_boxplot:
    num_iters = number_of_subjects
for sub_plot in range(2,1+num_iters):

    # Iterate through all the files starting with with and without haptic forces
    for group in range(len(hapticforces)):

        # The average of each amplitude signal is taken. num_0percent-num_50percent records the number of signals added together
        num_freq = np.zeros(4)
        A_Fmag = np.zeros((4,w_len))
        A_Fx = np.zeros((4,w_len))
        A_Fy = np.zeros((4,w_len))
        if split_signal:
            A_Fmag_highenergy = np.zeros((4,w_len))
            A_Fx_highenergy = np.zeros((4,w_len))
            A_Fy_highenergy = np.zeros((4,w_len))
            A_Fmag_lowenergy = np.zeros((4,w_len))
            A_Fx_lowenergy = np.zeros((4,w_len))
            A_Fy_lowenergy = np.zeros((4,w_len))
            if group == 1:
                break
        for subject_num in range(2,1+number_of_subjects): #iterate though subjects

            if (make_plot_each_sub==0) or (make_plot_each_sub==1 and sub_plot==subject_num):
                # Sets the correct subfile
                subname = "S0" + str(subject_num)
                subfile = "/OutputData_S" + str(subject_num)

                for freq in range(0,4):
                    if split_signal:
                        percent_height = 0.3
                        gravity = 9.81
                        mass = 1.0
                        max_energy = mass * gravity * percent_height * radius_options[freq]
                        mid_energy = mass * gravity * percent_height * radius_options[freq] / 2.0
                    for k in range(1,50): #iterate through trials
                        try:

                            # Open the trial files
                            trialfile = DIR+subname+subfile+"_Trial"+str(k)+"_Freq"+str(freq)+"_SL0_F"+str(group)+".csv"
                            # print(trialfile)
                            data = genfromtxt(trialfile,delimiter=',',dtype=float)
                            data = np.delete(data,0,0) # Deletes the first column of column names
                            print(trialfile)


                            # Get amplitudes for force and position signals
                            # Combine the x and y direction forces
                            y = np.sqrt(np.square(data[:,8])+np.square(data[:,9]))
                            A_Fmag_i,frq = calculate_amplitude(w,y,Fs,False,True)

                            # Use signal directly 1-x position 2- y position 8- x force 9- y force
                            y_Fx = data[:,8]
                            A_Fx_i,frq = calculate_amplitude(w,y_Fx,Fs,False,True)
                            y_Fy = data[:,9]
                            A_Fy_i,frq = calculate_amplitude(w,y_Fy,Fs,False,True)

                            A_Fmag[freq,:] += A_Fmag_i
                            A_Fx[freq,:] += A_Fx_i
                            A_Fy[freq,:] += A_Fy_i
                            num_freq[freq] += 1

                            if split_signal:
                                y_highenergy = []
                                y_Fx_highenergy = []
                                y_Fy_highenergy = []
                                y_lowenergy = []
                                y_Fx_lowenergy = []
                                y_Fy_lowenergy = []
                                count = 0
                                for row in range(0,data.shape[0]):
                                    if data[row,16]>max_energy*3:
                                        count += 1
                                    if data[row,16]>mid_energy:# and data[row,16]<max_energy*3:
                                        y_highenergy.append(np.sqrt(np.square(data[row,8])+np.square(data[row,9])))
                                        y_Fx_highenergy.append(data[row,8])
                                        y_Fy_highenergy.append(data[row,9])
                                    elif data[row,16]<mid_energy:
                                        y_lowenergy.append(np.sqrt(np.square(data[row,8])+np.square(data[row,9])))
                                        y_Fx_lowenergy.append(data[row,8])
                                        y_Fy_lowenergy.append(data[row,9])
                                print(count)
                                A_Fmag_highenergy_i,frq = calculate_amplitude(w,y_highenergy,Fs,False,True)
                                A_Fx_highenergy_i,frq = calculate_amplitude(w,y_Fx_highenergy,Fs,False,True)
                                A_Fy_highenergy_i,frq = calculate_amplitude(w,y_Fy_highenergy,Fs,False,True)

                                A_Fmag_lowenergy_i,frq = calculate_amplitude(w,y_lowenergy,Fs,False,True)
                                A_Fx_lowenergy_i,frq = calculate_amplitude(w,y_Fx_lowenergy,Fs,False,True)
                                A_Fy_lowenergy_i,frq = calculate_amplitude(w,y_Fy_lowenergy,Fs,False,True)

                                if np.min(A_Fmag_highenergy_i)>0:
                                    A_Fmag_highenergy[freq,:] += A_Fmag_highenergy_i
                                    A_Fx_highenergy[freq,:] += A_Fx_highenergy_i
                                    A_Fy_highenergy[freq,:] += A_Fy_highenergy_i

                                if np.min(A_Fmag_lowenergy_i)>0:
                                    A_Fmag_lowenergy[freq,:] += A_Fmag_lowenergy_i
                                    A_Fx_lowenergy[freq,:] += A_Fx_lowenergy_i
                                    A_Fy_lowenergy[freq,:] += A_Fy_lowenergy_i

                                # num_freq[freq] += 1

                        except:
                            pass

        if make_plots==1:
            # Take average of each signal
            for freq in range(0,4):
                if num_freq[freq]!=0:
                    dw = w[1]-w[0]
                    if split_signal:
                        A_Fmag_highenergy[freq,:] /= num_freq[freq]
                        A_Fx_highenergy[freq,:] /= num_freq[freq]
                        A_Fy_highenergy[freq,:] /= num_freq[freq]
                        E_signal = np.sum(np.square(A_Fmag_highenergy[freq,:]))*dw
                        A_Fmag_highenergy[freq,:] *= np.sqrt(1/E_signal)
                        E_signal = np.sum(np.square(A_Fx_highenergy[freq,:]))*dw
                        A_Fx_highenergy[freq,:] *= np.sqrt(1/E_signal)
                        E_signal = np.sum(np.square(A_Fy_highenergy[freq,:]))*dw
                        A_Fy_highenergy[freq,:] *= np.sqrt(1/E_signal)

                        A_Fmag_lowenergy[freq,:] /= num_freq[freq]
                        A_Fx_lowenergy[freq,:] /= num_freq[freq]
                        A_Fy_lowenergy[freq,:] /= num_freq[freq]
                        E_signal = np.sum(np.square(A_Fmag_lowenergy[freq,:]))*dw
                        A_Fmag_lowenergy[freq,:] *= np.sqrt(1/E_signal)
                        E_signal = np.sum(np.square(A_Fx_lowenergy[freq,:]))*dw
                        A_Fx_lowenergy[freq,:] *= np.sqrt(1/E_signal)
                        E_signal = np.sum(np.square(A_Fy_lowenergy[freq,:]))*dw
                        A_Fy_lowenergy[freq,:] *= np.sqrt(1/E_signal)

                    else:
                        A_Fmag[freq,:] /= num_freq[freq]
                        A_Fx[freq,:] /= num_freq[freq]
                        A_Fy[freq,:] /= num_freq[freq]

                        E_signal = np.sum(np.square(A_Fmag[freq,:]))*dw
                        A_Fmag[freq,:] *= np.sqrt(1/E_signal)

                        E_signal = np.sum(np.square(A_Fx[freq,:]))*dw
                        A_Fx[freq,:] *= np.sqrt(1/E_signal)

                        E_signal = np.sum(np.square(A_Fy[freq,:]))*dw
                        A_Fy[freq,:] *= np.sqrt(1/E_signal)

                        mag_list.append(A_Fmag[freq,:])
                        xy_xlist.append(A_Fx[freq,:])
                        xy_ylist.append(A_Fy[freq,:])


                    if freq_range_boxplot and group==0:
                        dw = w[1]-w[0]
                        for freq2 in range(0,4):
                            freq_list = []
                            for w_i in range(0,len(w)):
                                if (w[w_i] < freq_pendulum[freq2]+window) and (w[w_i] > freq_pendulum[freq2]-window):
                                    freq_list.append(A_Fx[freq,w_i])
                            energy_mat_fx[freq,freq2,sub_plot-2] += np.sum(np.square(freq_list))*dw

                            freq_list = []
                            for w_i in range(0,len(w)):
                                if (w[w_i] < freq_pendulum[freq2]+window) and (w[w_i] > freq_pendulum[freq2]-window):
                                    freq_list.append(A_Fy[freq,w_i])
                            energy_mat_fy[freq,freq2,sub_plot-2] += np.sum(np.square(freq_list))*dw


            if split_signal:# and group==1:
                for freq in range(0,4):
                    mag_list.append(A_Fmag_lowenergy[freq,:])
                    xy_xlist.append(A_Fx_lowenergy[freq,:])
                    xy_ylist.append(A_Fy_lowenergy[freq,:])
                for freq in range(0,4):
                    mag_list.append(A_Fmag_highenergy[freq,:])
                    xy_xlist.append(A_Fx_highenergy[freq,:])
                    xy_ylist.append(A_Fy_highenergy[freq,:])
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
        legend = ["Ball's resonant\nfrequency",
            'Low energy ball '+frequencylabels[0],
            'Low energy ball '+frequencylabels[1],
            'Low energy ball '+frequencylabels[2],
            'Low energy ball '+frequencylabels[3],
            'High energy ball '+frequencylabels[0],
            'High energy ball '+frequencylabels[1],
            'High energy ball '+frequencylabels[2],
            'High energy ball '+frequencylabels[3]]
        savename = 'xy_freq_split.png'
    else:
        legend = ["Ball's resonant\nfrequency",
            'forces '+hapticforces[0]+' '+frequencylabels[0],
            'forces '+hapticforces[0]+' '+frequencylabels[1],
            'forces '+hapticforces[0]+' '+frequencylabels[2],
            'forces '+hapticforces[0]+' '+frequencylabels[3],
            'forces '+hapticforces[1]+' '+frequencylabels[0],
            'forces '+hapticforces[1]+' '+frequencylabels[1],
            'forces '+hapticforces[1]+' '+frequencylabels[2],
            'forces '+hapticforces[1]+' '+frequencylabels[3]]
        savename = 'xy_freq.png'

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
    ####################################################################
    ######   Frequency magnitude spectrum plot for each frequency  #####
    ####################################################################
    ####################################################################
    for freq in range(0,4):
        title = 'Frequency Spectrum Plot for ' + frequencylabels[freq]
        if split_signal:
            savename = 'xy_freq_split_'+frequencylabels[freq]+'.png'
            legend = ["Ball's resonant\nfrequency",
                        'Low ball energy',
                        'High ball energy']
        else:
            savename = 'xy_freq_'+frequencylabels[0]+'.png'
            legend = ["Ball's resonant\nfrequency",
                        'forces '+hapticforces[0],
                        'forces '+hapticforces[1]]
        # colors_fq = [colors[freq],colors[freq+4]]
        colors_fq = [colors[0],colors[1]]
        # linestyles_fq = [linestyles[freq],linestyles[freq+4]]
        linestyles_fq = [linestyles[4],linestyles[5]]
        xdata = []
        ydata = []
        xdata.append(xy_xlist[freq])
        xdata.append(xy_xlist[freq+4])
        ydata.append(xy_ylist[freq])
        ydata.append(xy_ylist[freq+4])
        [fig,ax] = xy_spectrum(w,xdata,ydata,[freq_pendulum[freq]],title,xlabel,ylabel,legend,linestyles_fq,colors_fq,ymin,ymax,ymax_pend)
        fig.savefig(savename)

if make_plots and freq_range_boxplot:
    # energy_mat_fx = np.divide(energy_mat_fx,energy_num)
    data = []
    for freq2 in range(0,4):
        for freq in range(0,4):

            data.append(np.array(energy_mat_fx[freq,freq2,:]))

    # Plot parameters
    labels = ('0.5-0.5','0.5-1','0.5-1.5','0.5-2','1-0.5','1-1','1-1.5','1-2','1.5-0.5','1.5-1','1.5-1.5','1.5-2','2-0.5','2-1','2-1.5','2-2')
    # box_colors = ['#601A4A','#601A4A','#601A4A','#601A4A',
    #     '#EE442F','#EE442F','#EE442F','#EE442F',
    #     '#63ACBE','#63ACBE','#63ACBE','#63ACBE',
    #     '#006400','#006400','#006400','#006400']
    # box_alpha = [0.25,.5,.75,1,0.25,.5,.75,1,0.25,.5,.75,1,0.25,.5,.75,1] # sets how transparent to make the color

    box_colors = ['#601A4A','#EE442F','#63ACBE','#006400','#601A4A','#EE442F','#63ACBE','#006400','#601A4A','#EE442F','#63ACBE','#006400','#601A4A','#EE442F','#63ACBE','#006400']
    box_alpha = [0.25,0.25,0.25,0.25,.5,.5,.5,.5,.75,.75,.75,.75,1,1,1,1]

    figure_size = (6,3.55) # sets the size of the figure in inches
    xlabel = 'Frequency (Range (Hz) - Ball (Hz))'
    ylabel = 'Energy'

    title = 'X-dir Energy Around Different Frequencies-- range +/-' + str(window)
    [fig,ax]=make_boxplot(data,title,xlabel,ylabel,labels,box_colors,box_alpha,figure_size)
    fig.savefig('energy_fx.png')

    # energy_mat_fy = np.divide(energy_mat_fy,energy_num)
    data = []
    for freq2 in range(0,4):
        for freq in range(0,4):
            data.append(np.array(energy_mat_fy[freq,freq2,:]))

    title = 'Y-dir Energy Around Different Frequencies-- range +/-' + str(window)
    [fig2,ax2]=make_boxplot(data,title,xlabel,ylabel,labels,box_colors,box_alpha,figure_size)
    fig2.savefig('energy_fy.png')
plt.show()
