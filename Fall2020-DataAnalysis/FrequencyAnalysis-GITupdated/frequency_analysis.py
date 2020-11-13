import numpy as np
import scipy as sp
from numpy import fft, genfromtxt
from scipy import fft, arange, stats
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import csv

# This program analyses the frequency spectrum for each trial/participant.
# It saves and shows figures for the BioRob paper and saves the values of the peaks
# around the pendulum's resonant frequency into csv files
# "freq-metrics.csv", "freq-metrics-controls.csv", "freq-metrics-stroke.csv"
# for statistical analyses in R.
# This code expects the original data to be formatted and named properly

# Edit these variables before running
save_values = 1 # 0-do not save 1-save
make_plots = 1 # 0-do not make plots 1-make plots
DIR = "Z:" #set directory where data is mounted Ola- "/media/ola/Elements/R01prelim" Milli -"Z:"
number_of_subjects = 6 # in each group, controls and stroke
DT = 0.01
Fs = 1/DT
freq_pendulum = (1/(2*np.pi))*np.sqrt(9.81/0.07) #resonant frequency of pendulum

def calculate_amplitude(w,y):
    """
    Calculates the discrete time fourier transform and frequency amplitudes for the signal,
    places the values in frequency "bins" specified by vector w,
    and normalizes so the energy of the spectrum is 1
    Inputs:
        w - frequency values for "bins", note that it must start at 0
        y - discrete values of function
    Outputs:
        A - vector of signal amplitude at each frequency in w (normalized)
    """

    # Use FFT to calculate the frequency spectrum for the given data length
    n = len(y) # length of the signal
    # Y = abs(np.fft.fft(y, n=n, norm='ortho')) #Normalize by signal length
    Y = abs(np.fft.fft(y, n=n)) #No normalization
    frq = np.fft.fftfreq(n, d=1/Fs) # Determine the corresponding frequency values for the fft output
    # Choose the one sided frequencies of interest
    frq = frq[range(n//2)] # one side frequency range
    Y = Y[range(n//2)]

    # Place values into prespecified frequency bins, each with range (w[w_index]-(dw/2),w[w_index]+(dw/2)]
    w_len = len(w)
    dw = w[1]-w[0]
    w_index = 0
    A = np.zeros(w_len)
    num = np.zeros(w_len)
    Y_len = len(Y)
    for i in range(0,Y_len): # iterate through the frequency outputs of fft
        # increase the index of the frequency if the frequency is greater than the frequency corresponding to the index of the current "bin"
        if frq[i]>(w[w_index]+(dw/2)):
            w_index += 1
        # print('i', i, 'frq[i]', frq[i], 'w_index', w_index, '(w[w_index]+(dw/2))', (w[w_index]+(dw/2)) )
        # print('w_index', w_index, 'w_len', w_len)

        # Stop once the index is too large. This is particilarly important if the bins do not go up the nyquist frequency
        if (w_index>=w_len):
            break
        # Add the amplitude to the "bin"
        A[w_index] += Y[i]
        num[w_index] += 1

    # Take the average of each "bin"
    for i in range(w_len):
        if (num[i]>0):
            A[i] = A[i]/num[i]

    # Calculate the energy of the amplitude signal
    E_signal = np.sum(np.square(A[0:w_len]))*dw
    # print('Signal energy: ',E_signal)

    # Calculate new signal normalized to have an energy of 1
    A *= np.sqrt(1/E_signal)
    E_signal = np.sum(np.square(A[0:w_len]))*dw
    # print('Signal energy: ',E_signal) # should equal to 1

    # # Plots transformation
    # t = arange(0,20,DT) # time vector
    # plt.subplot(3,1,1)
    # plt.plot(t,y)
    # plt.xlabel('Time')
    # plt.ylabel('Amplitude')
    # plt.subplot(3,1,2)
    # plt.plot(frq,Y,'r') # plotting the spectrum
    # plt.xlabel('Freq (Hz)')
    # plt.ylabel('Amplitude')
    # plt.subplot(3,1,3)
    # plt.plot(w,A,'b') # plotting the spectrum
    # plt.xlabel('Freq (Hz)')
    # plt.ylabel('Binned amplitude')
    # plt.show()
    return A


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
supportlevels = ['0%Max','20%Max','50%Max']
tasks = ['Task1','Task2','Task3','Task4','Task5']
subjects = ['Subject1','Subject2','Subject3','Subject4','Subject5','Subject6','Subject7','Subject8','Subject9','Subject10','Subject11','Subject12','Subject13','Subject14','Subject15','Subject16','Subject17','Subject18']
colors = ['#601A4A','#EE442F','#63ACBE',"#E0E0E0"] # for plots
group_name = ['Controls','Stroke']
linestyles = ['--','-']
markerstyles = ['o','D']

# Store data for statistical analysis
if save_values==1:
    file = DIR+"freq-metrics.csv"
    columns = ["Impairment","Subject","SupportLevel","Task","Trial","FmagPeak","DiffPeak","Ratio","FxPeak","FyPeak","PxPeak","PyPeak","RatioFx","RatioFy","FxHigh"]
    with open(file,'w') as csvfile:
        testwriter = csv.writer(csvfile,delimiter=',')
        testwriter.writerow(columns)
    file_controls = DIR+"freq-metrics-controls.csv"
    columns = ["Subject","SupportLevel","Task","Trial","FmagPeak","DiffPeak","Ratio","FxPeak","FyPeak","PxPeak","PyPeak","RatioFx","RatioFy","FxHigh"]
    with open(file_controls,'w') as csvfile:
        testwriter = csv.writer(csvfile,delimiter=',')
        testwriter.writerow(columns)
    file_stroke = DIR+"freq-metrics-stroke.csv"
    with open(file_stroke,'w') as csvfile:
        testwriter = csv.writer(csvfile,delimiter=',')
        testwriter.writerow(columns)

# Initialize arrays for storing row = support level col = subject
low_freq_Fx = np.zeros((3,number_of_subjects))
high_freq_Fx = np.zeros((3,number_of_subjects))
low_freq_Fy = np.zeros((3,number_of_subjects))
high_freq_Fy = np.zeros((3,number_of_subjects))

# Initialize arrays for creating boxplot with
# row = subject
# col = in order [stroke0%; control0%; stroke30%; control30%; stroke50%; control50%]
data_boxplot = np.zeros((number_of_subjects,6))
data_boxplot_low = np.zeros((number_of_subjects,6))

# Iterate through all the files
for group in range(len(group_name)):
    # The average of each amplitude signal is taken. num_0percent-num_50percent records the number of signals added together
    num_0percent = 0
    num_20percent = 0
    num_50percent = 0

    for subject_num in range(1,1+number_of_subjects): #iterate though subjects

        # Sets the correct subfile depending on which group the subject is in
        if (group_name[group] == 'Stroke'):
            subname = "S0" + str(subject_num)
        else:
            subject_num += 10
            subname = "S" + str(subject_num)
        subfile = "/OutputData_S" + str(subject_num)

        # For plotting individual subject figures
        # plt.figure(subject_num)
        # plt.title('Subject '+str(subject_num)+' Frequency Spectrum')
        for j in range(0,3): #iterate through support levels
            # num_supportlevel = 0 # keeps track of the number of subjects at that support level for individual subject analyses
            for k in range(1,46): #iterate through trials
                for i in range(1,6): #iterate through tasks
                    try:
                        # Open the trial files
                        trialfile = DIR+subname+subfile+"_Task"+str(i)+"_SL"+str(j)+"_Trial"+ str(k)+".csv"
                        data = genfromtxt(trialfile,delimiter=',',dtype=float)
                        data = np.delete(data,0,0) # Deletes the first column of column names
                        print(trialfile)
                        # Determine the indecies for the start and end of trial
                        delay = data[-1,0]-20.0 # Gets the final time to calculate the delay in data collection
                        delayInd = int(delay/DT)
                        score = data[-1,7]
                        # gets end of trial by moving backwards until the score is less than 20
                        trial_complete_index = int(data.shape[0]-1)
                        if score==20:
                            while data[trial_complete_index,7]==20:
                                trial_complete_index-=1
                        # print('index: ',trial_complete_index, 'time: ', data[trial_complete_index,0]-delay)

                        # Get amplitudes for force and position signals
                        # Combine the x and y direction forces
                        y = np.sqrt(np.square(data[delayInd:trial_complete_index,8])+np.square(data[delayInd:trial_complete_index,9]))
                        A_Fmag = calculate_amplitude(w,y)
                        # Use signal directly 1-x position 2- y position 8- x force 9- y force
                        y_Fx = data[delayInd:trial_complete_index,8]
                        A_Fx = calculate_amplitude(w,y_Fx)
                        y_Fy = data[delayInd:trial_complete_index,9]
                        A_Fy = calculate_amplitude(w,y_Fy)
                        x_acc = data[delayInd:trial_complete_index,12]
                        y_Px = y_Fx + x_acc
                        # y = data[delayInd:trial_complete_index,1] position
                        A_Px = calculate_amplitude(w,y_Px)
                        y_acc = data[delayInd:trial_complete_index,12]
                        y_Py = y_Fy + y_acc
                        # y = data[delayInd:trial_complete_index,2] position
                        A_Py = calculate_amplitude(w,y_Py)
                        y_Fmag_net = np.sqrt(np.square(y_Px)+np.square(y_Py))
                        A_Fmag_net = calculate_amplitude(w,y_Fmag_net)
                        A_Dx = A_Fx - A_Px
                        A_Dy = A_Fy - A_Py
                        A_Fmag_diff = A_Fmag - A_Fmag_net

                        if make_plots==1:
                            # For combining participants for plots
                            if (j==0):
                                if (num_0percent==0):
                                    A_Fmag_0percent = A_Fmag
                                    A_Fx_0percent = A_Fx
                                    A_Fy_0percent = A_Fy
                                    A_Px_0percent = A_Px
                                    A_Py_0percent = A_Py
                                    A_Dx_0percent = A_Dx
                                    A_Dy_0percent = A_Dy
                                    A_Fmag_net_0percent = A_Fmag_net
                                    A_Fmag_diff_0percent = A_Fmag_diff
                                else:
                                    A_Fmag_0percent += A_Fmag
                                    A_Fx_0percent += A_Fx
                                    A_Fy_0percent += A_Fy
                                    A_Px_0percent += A_Px
                                    A_Py_0percent += A_Py
                                    A_Dx_0percent += A_Dx
                                    A_Dy_0percent += A_Dy
                                    A_Fmag_net_0percent += A_Fmag_net
                                    A_Fmag_diff_0percent += A_Fmag_diff
                                num_0percent += 1
                            if (j==1):
                                if (num_20percent==0):
                                    A_Fmag_20percent = A_Fmag
                                    A_Fx_20percent = A_Fx
                                    A_Fy_20percent = A_Fy
                                    A_Px_20percent = A_Px
                                    A_Py_20percent = A_Py
                                    A_Dx_20percent = A_Dx
                                    A_Dy_20percent = A_Dy
                                    A_Fmag_net_20percent = A_Fmag_net
                                    A_Fmag_diff_20percent = A_Fmag_diff
                                else:
                                    A_Fmag_20percent += A_Fmag
                                    A_Fx_20percent += A_Fx
                                    A_Fy_20percent += A_Fy
                                    A_Px_20percent += A_Px
                                    A_Py_20percent += A_Py
                                    A_Dx_20percent += A_Dx
                                    A_Dy_20percent += A_Dy
                                    A_Fmag_net_20percent += A_Fmag_net
                                    A_Fmag_diff_20percent += A_Fmag_diff
                                num_20percent += 1
                            if (j==2):
                                if (num_50percent==0):
                                    A_Fmag_50percent = A_Fmag
                                    A_Fx_50percent = A_Fx
                                    A_Fy_50percent = A_Fy
                                    A_Px_50percent = A_Px
                                    A_Py_50percent = A_Py
                                    A_Dx_50percent = A_Dx
                                    A_Dy_50percent = A_Dy
                                    A_Fmag_net_50percent = A_Fmag_net
                                    A_Fmag_diff_50percent = A_Fmag_diff
                                else:
                                    A_Fmag_50percent += A_Fmag
                                    A_Fx_50percent += A_Fx
                                    A_Fy_50percent += A_Fy
                                    A_Px_50percent += A_Px
                                    A_Py_50percent += A_Py
                                    A_Dx_50percent += A_Dx
                                    A_Dy_50percent += A_Dy
                                    A_Fmag_net_50percent += A_Fmag_net
                                    A_Fmag_diff_50percent += A_Fmag_diff
                                num_50percent += 1

                        if save_values==1 or make_plots==1:
                            # Write max value of peak to csv file for statistical analysis
                            resonant_peak_Fmag = 0
                            resonant_peak_Fmag_diff = 0
                            resonant_peak_Fx = 0
                            resonant_peak_Fy = 0
                            resonant_peak_Px = 0
                            resonant_peak_Py = 0
                            low_freq_energy = 0
                            high_freq_energy = 0
                            low_freq_energy_Fx = 0
                            high_freq_energy_Fx = 0
                            low_freq_energy_Fy = 0
                            high_freq_energy_Fy = 0
                            for w_i in range(0,len(w)):
                                # Get peak for force magnitude around 2*freq_pendulum
                                if (w[w_i] < 2*freq_pendulum+1) and (w[w_i] > 2*freq_pendulum-1):
                                    if (A_Fmag[w_i]>resonant_peak_Fmag):
                                        resonant_peak_Fmag = A_Fmag[w_i]
                                    if (A_Fmag_diff[w_i]>resonant_peak_Fmag_diff):
                                        resonant_peak_Fmag_diff = A_Fmag_diff[w_i]
                                # Get peak around force x and y
                                if (w[w_i] < freq_pendulum+1) and (w[w_i] > freq_pendulum-1):
                                    if (A_Fx[w_i]>resonant_peak_Fx):
                                        resonant_peak_Fx = A_Fx[w_i]
                                    if (A_Fy[w_i]>resonant_peak_Fy):
                                        resonant_peak_Fy = A_Fy[w_i]
                                    if (A_Px[w_i]>resonant_peak_Px):
                                        resonant_peak_Px = A_Px[w_i]
                                    if (A_Py[w_i]>resonant_peak_Py):
                                        resonant_peak_Py = A_Py[w_i]
                                if (w[w_i]<1):
                                    low_freq_energy_Fx += A_Fx[w_i]
                                    low_freq_energy_Fy += A_Fy[w_i]
                                    low_freq_energy += A_Fmag[w_i]
                                    last_under_1_Hz = w_i
                                else:
                                    high_freq_energy_Fx += A_Fx[w_i]
                                    high_freq_energy_Fy += A_Fy[w_i]
                                    high_freq_energy += A_Fmag[w_i]

                            # Find ratio of the energy of the signal at high frequencies vs. low
                            # no need to multiply by dw because it is a ratio and dw would cancel out
                            ratio_energy = high_freq_energy/low_freq_energy
                            ratio_energy_Fx = high_freq_energy_Fx/low_freq_energy_Fx
                            ratio_energy_Fy = high_freq_energy_Fy/low_freq_energy_Fy
                            # store for ratio plot
                            # ratio_index used to that subjects 1-6 and 11-16 each correspond to indecies 0-5
                            ratio_index = subject_num-1
                            supportlevel_index = j*2
                            if (group_name[group] == 'Controls'):
                                ratio_index -= 10
                                supportlevel_index = 1 + j*2
                            low_freq_Fx[j,ratio_index] += low_freq_energy_Fx  # support level is the row and task is the column
                            high_freq_Fx[j,ratio_index] += high_freq_energy_Fx  # support level is the row and task is the column
                            low_freq_Fy[j,ratio_index] += low_freq_energy_Fy  # support level is the row and task is the column
                            high_freq_Fy[j,ratio_index] += high_freq_energy_Fy  # support level is the row and task is the column

                            dw = w[1]-w[0]
                            Fx_high_content = np.sum(np.square(A_Fx[0:last_under_1_Hz]))*dw
                            data_boxplot[ratio_index,supportlevel_index] += np.sum(np.square(A_Fx[last_under_1_Hz+1:]))*dw
                            data_boxplot_low[ratio_index,supportlevel_index] += np.sum(np.square(A_Fx[0:last_under_1_Hz]))*dw
                            # print(data_boxplot_low[ratio_index,supportlevel_index],data_boxplot[ratio_index,supportlevel_index],data_boxplot[ratio_index,supportlevel_index]+data_boxplot_low[ratio_index,supportlevel_index])

                        if save_values==1:
                            row = [resonant_peak_Fmag, resonant_peak_Fmag_diff, ratio_energy, resonant_peak_Fx, resonant_peak_Fy, resonant_peak_Px, resonant_peak_Py, ratio_energy_Fx, ratio_energy_Fy,Fx_high_content]
                            # inserts the trial/subject information before metrics
                            row.insert(0,k) # trial number
                            row.insert(0,tasks[i-1])
                            row.insert(0,supportlevels[j])
                            row.insert(0,subjects[subject_num-1])
                            if (group_name[group] == 'Stroke'):
                                with open (file_stroke,'a') as csvfile:
                                    testwriter_stroke = csv.writer(csvfile,delimiter=',')
                                    testwriter_stroke.writerow(row)
                            else:
                                with open (file_controls,'a') as csvfile:
                                    testwriter_controls = csv.writer(csvfile,delimiter=',')
                                    testwriter_controls.writerow(row)
                            row.insert(0,group_name[group])
                            with open (file,'a') as csvfile:
                                testwriter = csv.writer(csvfile,delimiter=',')
                                testwriter.writerow(row)

                        # # For making individual plot
                        # if (num_supportlevel==0):
                        #     A_support = A
                        # else:
                        #     A_support += A
                        # num_supportlevel += 1
                    except:
                        pass
                # For individual subject plot
                # A_support /= num_supportlevel
                # plt.plot(w[wi_0:w_len],A_support[wi_0:w_len],linestyle=linestyles[group],color=colors[j],label=supportlevels[j]) # plotting the spectrum (start at )
            # plt.legend(loc="upper right")
            # plt.xlabel('Freq (Hz)')
            # plt.ylabel('Amplitude')
            # plt.xscale('log')
            # plt.yscale('log')
            # plt.savefig('sub'+str(subject_num)+'_freq.png')
    if make_plots==1:
        # Take average of each signal
        A_Fmag_0percent /= num_0percent
        A_Fmag_20percent /= num_20percent
        A_Fmag_50percent /= num_50percent
        A_Fx_0percent /= num_0percent
        A_Fx_20percent /= num_20percent
        A_Fx_50percent /= num_50percent
        A_Fy_0percent /= num_0percent
        A_Fy_20percent /= num_20percent
        A_Fy_50percent /= num_50percent
        A_Px_0percent /= num_0percent
        A_Px_20percent /= num_20percent
        A_Px_50percent /= num_50percent
        A_Py_0percent /= num_0percent
        A_Py_20percent /= num_20percent
        A_Py_50percent /= num_50percent
        A_Dx_0percent /= num_0percent
        A_Dx_20percent /= num_20percent
        A_Dx_50percent /= num_50percent
        A_Dy_0percent /= num_0percent
        A_Dy_20percent /= num_20percent
        A_Dy_50percent /= num_50percent
        A_Fmag_net_0percent /= num_0percent
        A_Fmag_net_20percent /= num_20percent
        A_Fmag_net_50percent /= num_50percent
        A_Fmag_diff_0percent /= num_0percent
        A_Fmag_diff_20percent /= num_20percent
        A_Fmag_diff_50percent /= num_50percent
        # Divide each value by 15 because there are 15 instances of each supporlevel/subject combination
        num_per_support = 3*5
        low_freq_Fx /= num_per_support
        high_freq_Fx /= num_per_support
        low_freq_Fy /= num_per_support
        high_freq_Fy /= num_per_support

        # Figures with both stroke and controls

        # # frequency ratio plot
        # figure_size = (8,3.55) # sets the size of the figure in inches
        # plt.figure(25, figsize=figure_size, dpi=150) #,figsize=(6.0, 4.0), dpi=80)
        # if (group==0):
        #     # create subplot
        #     fig_ratio, ax_ratio =plt.subplots(nrows=1, ncols=2, sharey='row', squeeze=True, figsize=figure_size, dpi=150)
        #     fig_ratio.subplots_adjust(hspace=0.05)
        #     fig_ratio.subplots_adjust(wspace=0.05)
        #     # ymin = 0.0003
        #     # ymax = 5
        #
        #     plot_num = 0
        #     ax_ratio[plot_num].grid(True, color=colors[3])
        #     l1 = ax_ratio[plot_num].scatter(low_freq_Fx[0,:],high_freq_Fx[0,:],marker=markerstyles[group],color=colors[0])
        #     l2 = ax_ratio[plot_num].scatter(low_freq_Fx[1,:],high_freq_Fx[1,:],marker=markerstyles[group],color=colors[1])
        #     l3 = ax_ratio[plot_num].scatter(low_freq_Fx[2,:],high_freq_Fx[2,:],marker=markerstyles[group],color=colors[2])
        #     # set plot parameters
        #     # ax_force[plot_num].set_xlim((10^-1,w[w_len-1]))
        #     # ax_force[plot_num].set_ylim(ymin,ymax)
        #     ax_ratio[plot_num].text(.15,.93,'X-Component',horizontalalignment='left',transform=ax_ratio[plot_num].transAxes, fontsize=10)
        #     ax_ratio[plot_num].text(.05,.93,'A.',horizontalalignment='left',transform=ax_ratio[plot_num].transAxes, fontsize=10, fontweight='bold')
        #
        #     plot_num = 1
        #     ax_ratio[plot_num].grid(True, color=colors[3])
        #     ax_ratio[plot_num].scatter(low_freq_Fy[0,:],high_freq_Fy[0,:],marker=markerstyles[group],color=colors[0])
        #     ax_ratio[plot_num].scatter(low_freq_Fy[1,:],high_freq_Fy[1,:],marker=markerstyles[group],color=colors[1])
        #     ax_ratio[plot_num].scatter(low_freq_Fy[2,:],high_freq_Fy[2,:],marker=markerstyles[group],color=colors[2])
        #     # # set plot parameters
        #     # # ax_force[plot_num].set_xlim((10^-1,w[w_len-1]))
        #     # # ax_force[plot_num].set_ylim(ymin,ymax)
        #     ax_ratio[plot_num].text(.15,.93,'Y-Component',horizontalalignment='left',transform=ax_ratio[plot_num].transAxes, fontsize=10)
        #     ax_ratio[plot_num].text(.05,.93,'B.',horizontalalignment='left',transform=ax_ratio[plot_num].transAxes, fontsize=10, fontweight='bold')
        #
        # if (group==1):
        #     plot_num = 0
        #     l4 = ax_ratio[plot_num].scatter(low_freq_Fx[0,:],high_freq_Fx[0,:],marker=markerstyles[group],color=colors[0])
        #     l5 = ax_ratio[plot_num].scatter(low_freq_Fx[1,:],high_freq_Fx[1,:],marker=markerstyles[group],color=colors[1])
        #     l6 = ax_ratio[plot_num].scatter(low_freq_Fx[2,:],high_freq_Fx[2,:],marker=markerstyles[group],color=colors[2])
        #     for label in (ax_force[plot_num].get_xticklabels() + ax_force[plot_num].get_yticklabels()):
        #         label.set_fontsize(8)
        #
        #     plot_num = 1
        #     ax_ratio[plot_num].scatter(low_freq_Fy[0,:],high_freq_Fy[0,:],marker=markerstyles[group],color=colors[0])
        #     ax_ratio[plot_num].scatter(low_freq_Fy[1,:],high_freq_Fy[1,:],marker=markerstyles[group],color=colors[1])
        #     ax_ratio[plot_num].scatter(low_freq_Fy[2,:],high_freq_Fy[2,:],marker=markerstyles[group],color=colors[2])
        #     for label in (ax_force[plot_num].get_xticklabels() + ax_force[plot_num].get_yticklabels()):
        #         label.set_fontsize(8)
        #
        #     # create titles
        #     fig_ratio.text(0.5, 0.91, 'High-Low Frequency Ratio', ha='center', fontsize=10, fontweight='bold')
        #     # fig_force.suptitle('Force Frequency Spectrums', fontsize=10, fontweight='bold')
        #     fig_ratio.text(0.5, 0.01, 'Low-frequency Energy', ha='center', fontsize=10)
        #     fig_ratio.text(0.065, 0.5, 'High-frequency Energy', va='center', rotation='vertical', fontsize=10)
        #     line_labels = [group_name[0]+' '+supportlevels[0],group_name[0]+' '+supportlevels[1],group_name[0]+' '+supportlevels[2],group_name[1]+' '+supportlevels[0],group_name[1]+' '+supportlevels[1],group_name[1]+' '+supportlevels[2]]
        #     fig_ratio.legend([l1, l2, l3, l4, l5, l6],labels=line_labels,loc="center right", fontsize=9)
        #     fig_ratio.subplots_adjust(right=0.78)
        #     fig_ratio.savefig('freq_ratio.pdf')
        #
        #
        # # frequency spectrum plot
        # wi_0 = 0 # starting index (usually either 0 or 1)
        # figure_size = (8,3.55) # inches
        # plt.figure(20, figsize=figure_size, dpi=150) #,figsize=(6.0, 4.0), dpi=80)
        # if (group==0):
        #     # create subplot
        #     fig_force, ax_force =plt.subplots(nrows=1, ncols=2, sharey='row', squeeze=True, figsize=figure_size, dpi=150)
        #     fig_force.subplots_adjust(hspace=0.05)
        #     fig_force.subplots_adjust(wspace=0.05)
        #     ymin = 0.0003
        #     ymax = 5
        #     ymax_pend = 2
        #     x1, x2, y1, y2 = 1.2, 2.75, 0.19, 1.75 # specify the limits for inset axis
        #     zoom = 2.5
        #
        #     plot_num = 0
        #     l0 = ax_force[plot_num].plot([freq_pendulum, freq_pendulum],[ymin, ymax_pend],linestyle=':',color='k')
        #     l1 = ax_force[plot_num].plot(w[wi_0:w_len],A_Fx_0percent[wi_0:w_len],linestyle=linestyles[group],color=colors[0],label=(group_name[group]+' '+supportlevels[0]))
        #     l2 = ax_force[plot_num].plot(w[wi_0:w_len],A_Fx_20percent[wi_0:w_len],linestyle=linestyles[group],color=colors[1],label=(group_name[group]+' '+supportlevels[1]))
        #     l3 = ax_force[plot_num].plot(w[wi_0:w_len],A_Fx_50percent[wi_0:w_len],linestyle=linestyles[group],color=colors[2],label=(group_name[group]+' '+supportlevels[2]))
        #     # set plot parameters
        #     ax_force[plot_num].set_xscale('log')
        #     ax_force[plot_num].set_yscale('log')
        #     ax_force[plot_num].set_xlim((10^-1,w[w_len-1]))
        #     ax_force[plot_num].set_ylim(ymin,ymax)
        #     ax_force[plot_num].grid(True, color=colors[3])
        #     ax_force[plot_num].text(.15,.93,'X-Component',horizontalalignment='left',transform=ax_force[plot_num].transAxes, fontsize=10)
        #     ax_force[plot_num].text(.05,.93,'A.',horizontalalignment='left',transform=ax_force[plot_num].transAxes, fontsize=10, fontweight='bold')
        #     # create inset axes
        #     axins_for0 = zoomed_inset_axes(ax_force[plot_num], zoom, loc='lower left')#,bbox_to_anchor=(0.5, 8, 4, 4)) # bbox_to_anchor=(10,8)) # zoom-factor: 2.5, location: upper-left
        #     axins_for0.set_xlim(x1, x2) # apply the x-limits
        #     axins_for0.set_ylim(y1, y2) # apply the y-limits
        #     plt.tick_params(
        #         axis='both',
        #         which='both',
        #         bottom=False,
        #         left=False,
        #         right=False,
        #         top=False,
        #         labelleft=False,
        #         labelbottom=False)
        #     axins_for0.set_xscale('log')
        #     axins_for0.set_yscale('log')
        #     axins_for0.grid(True, color=colors[3], which='both')
        #     mark_inset(ax_force[plot_num], axins_for0, loc1=2, loc2=4, fc="none", ec="0.5",edgecolor='#000000', facecolor='#000000', color='#000000')
        #     axins_for0.plot([freq_pendulum, freq_pendulum],[ymin, ymax],linestyle=':',color='k')
        #
        #     plot_num = 1
        #     ax_force[plot_num].plot([freq_pendulum, freq_pendulum],[ymin, ymax_pend],linestyle=':',color='k')
        #     ax_force[plot_num].plot(w[wi_0:w_len],A_Fy_0percent[wi_0:w_len],linestyle=linestyles[group],color=colors[0],label=(group_name[group]+' '+supportlevels[0]))
        #     ax_force[plot_num].plot(w[wi_0:w_len],A_Fy_20percent[wi_0:w_len],linestyle=linestyles[group],color=colors[1],label=(group_name[group]+' '+supportlevels[1]))
        #     ax_force[plot_num].plot(w[wi_0:w_len],A_Fy_50percent[wi_0:w_len],linestyle=linestyles[group],color=colors[2],label=(group_name[group]+' '+supportlevels[2]))
        #     # set plot parameters
        #     ax_force[plot_num].set_xscale('log')
        #     ax_force[plot_num].set_yscale('log')
        #     ax_force[plot_num].set_xlim((10^-1,w[w_len-1]))
        #     ax_force[plot_num].set_ylim(ymin,ymax)
        #     ax_force[plot_num].grid(True, color=colors[3])
        #     ax_force[plot_num].text(.15,.93,'Y-Component',horizontalalignment='left',transform=ax_force[plot_num].transAxes, fontsize=10)
        #     ax_force[plot_num].text(.05,.93,'B.',horizontalalignment='left',transform=ax_force[plot_num].transAxes, fontsize=10, fontweight='bold')
        #     # create inset axes
        #     axins_for1 = zoomed_inset_axes(ax_force[plot_num], zoom, loc='lower left')#,bbox_to_anchor=(0.5, 8, 4, 4)) # bbox_to_anchor=(10,8)) # zoom-factor: 2.5, location: upper-left
        #     axins_for1.set_xlim(x1, x2) # apply the x-limits
        #     axins_for1.set_ylim(y1, y2) # apply the y-limits
        #     plt.tick_params(
        #         axis='both',
        #         which='both',
        #         bottom=False,
        #         left=False,
        #         right=False,
        #         top=False,
        #         labelleft=False,
        #         labelbottom=False)
        #     axins_for1.set_xscale('log')
        #     axins_for1.set_yscale('log')
        #     axins_for1.grid(True, color=colors[3], which='both')
        #     mark_inset(ax_force[plot_num], axins_for1, loc1=2, loc2=4, fc="none", ec="0.5",edgecolor='#000000', facecolor='#000000', color='#000000')
        #     axins_for1.plot([freq_pendulum, freq_pendulum],[ymin, ymax],linestyle=':',color='k')
        #
        # axins_for0.plot(w[wi_0:w_len],A_Fx_0percent[wi_0:w_len],linestyle=linestyles[group],color=colors[0],label=(group_name[group]+' '+supportlevels[0]))
        # axins_for0.plot(w[wi_0:w_len],A_Fx_20percent[wi_0:w_len],linestyle=linestyles[group],color=colors[1],label=(group_name[group]+' '+supportlevels[1]))
        # axins_for0.plot(w[wi_0:w_len],A_Fx_50percent[wi_0:w_len],linestyle=linestyles[group],color=colors[2],label=(group_name[group]+' '+supportlevels[2]))
        # axins_for1.plot(w[wi_0:w_len],A_Fy_0percent[wi_0:w_len],linestyle=linestyles[group],color=colors[0],label=(group_name[group]+' '+supportlevels[0]))
        # axins_for1.plot(w[wi_0:w_len],A_Fy_20percent[wi_0:w_len],linestyle=linestyles[group],color=colors[1],label=(group_name[group]+' '+supportlevels[1]))
        # axins_for1.plot(w[wi_0:w_len],A_Fy_50percent[wi_0:w_len],linestyle=linestyles[group],color=colors[2],label=(group_name[group]+' '+supportlevels[2]))
        #
        # if (group==1):
        #     plot_num = 0
        #     l4 = ax_force[plot_num].plot(w[wi_0:w_len],A_Fx_0percent[wi_0:w_len],linestyle=linestyles[group],color=colors[0],label=(group_name[group]+' '+supportlevels[0]))
        #     l5 = ax_force[plot_num].plot(w[wi_0:w_len],A_Fx_20percent[wi_0:w_len],linestyle=linestyles[group],color=colors[1],label=(group_name[group]+' '+supportlevels[1]))
        #     l6 = ax_force[plot_num].plot(w[wi_0:w_len],A_Fx_50percent[wi_0:w_len],linestyle=linestyles[group],color=colors[2],label=(group_name[group]+' '+supportlevels[2]))
        #     for label in (ax_force[plot_num].get_xticklabels() + ax_force[plot_num].get_yticklabels()):
        #         label.set_fontsize(8)
        #
        #     plot_num = 1
        #     ax_force[plot_num].plot(w[wi_0:w_len],A_Fy_0percent[wi_0:w_len],linestyle=linestyles[group],color=colors[0],label=(group_name[group]+' '+supportlevels[0]))
        #     ax_force[plot_num].plot(w[wi_0:w_len],A_Fy_20percent[wi_0:w_len],linestyle=linestyles[group],color=colors[1],label=(group_name[group]+' '+supportlevels[1]))
        #     ax_force[plot_num].plot(w[wi_0:w_len],A_Fy_50percent[wi_0:w_len],linestyle=linestyles[group],color=colors[2],label=(group_name[group]+' '+supportlevels[2]))
        #     for label in (ax_force[plot_num].get_xticklabels() + ax_force[plot_num].get_yticklabels()):
        #         label.set_fontsize(8)
        #
        #     # create titles
        #     fig_force.text(0.5, 0.91, 'Force Frequency Spectrums', ha='center', fontsize=10, fontweight='bold')
        #     # fig_force.suptitle('Force Frequency Spectrums', fontsize=10, fontweight='bold')
        #     fig_force.text(0.5, 0.01, 'Frequency (Hz)', ha='center', fontsize=10)
        #     fig_force.text(0.065, 0.5, 'Normalized Frequency Amplitude', va='center', rotation='vertical', fontsize=10)
        #     # fig_force.legend(loc="upper right", fontsize=9)
        #     line_labels = ["Ball's resonant\nfrequency",group_name[0]+' '+supportlevels[0],group_name[0]+' '+supportlevels[1],group_name[0]+' '+supportlevels[2],group_name[1]+' '+supportlevels[0],group_name[1]+' '+supportlevels[1],group_name[1]+' '+supportlevels[2]]
        #     fig_force.legend([l0, l1, l2, l3, l4, l5, l6],labels=line_labels,loc="center right", fontsize=9)
        #     # fig.legend([l1, l2, l3, l4, l5, l6],labels=line_labels,loc="center right", fontsize=9)
        #     fig_force.subplots_adjust(right=0.78)
        #     # label_x = 8
        #     # label_y = 0.0007
        #     # arrow_x = freq_pendulum*2
        #     # arrow_y = 0.003
        #     #
        #     # label_x = 18
        #     # label_y = 0.1
        #     # arrow_x = freq_pendulum*2
        #     # arrow_y = 0.6
        #     # arrow_properties = dict(
        #     #     facecolor="black", width=0.5,
        #     #     headwidth=4, shrink=0.1)
        #     # ax_force.annotate(
        #     #     "Pendulum's\nresonant\nfrequency", xy=(arrow_x, arrow_y),
        #     #     xytext=(label_x, label_y),
        #     #     arrowprops=arrow_properties, ha='center')
        #     fig_force.savefig('freq_forces.pdf')
if make_plots==1:

    # frequency ratio box plot
    num_per_support = 3*5
    data_boxplot /= num_per_support

    figure_size = (8,3.55) # sets the size of the figure in inches

    fig, ax1 = plt.subplots(figsize=figure_size)
    print(data_boxplot)
    labels = ['Stroke-0%','Controls-0%','Stroke-20%','Controls-20%','Stroke-50%','Controls-50%']
    medianprops = dict(linewidth=2.5, color='black')
    bp = ax1.boxplot(data_boxplot, notch=0, medianprops=medianprops, labels=labels)#, patch_artist=True,boxprops=dict(facecolor=color_combine, color=c))
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='+')
    # Add a horizontal grid to the plot, but make it very light in color
    # so we can use it for reading data values but not be distracting
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                   alpha=0.5)

    # Hide these grid behind plot objects
    ax1.set_axisbelow(True)
    ax1.set_title('High-Frequency Content For Each Loading Level')#, fontsize=10, fontweight='bold')
    ax1.set_xlabel('Loading Level (% Max)')
    ax1.set_ylabel('High-Frequency Content')
    # for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
    #     label.set_fontsize(8)

    box_colors = ['#601A4A', '#601A4A','#EE442F','#EE442F', '#63ACBE', '#63ACBE']
    box_alpha = [None,0.3,None,0.3,None,0.3]
    medians = np.empty(6)
    for i in range(6):
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        box_coords = np.column_stack([boxX, boxY])
        ax1.add_patch(Polygon(box_coords, facecolor=box_colors[i], alpha=box_alpha[i]))

    x_coordinates = np.array([1,2,3,4,5,6])
    y_coordinates = np.mean(data_boxplot,axis=0)
    # print(np.mean(data_boxplot_low,axis=0))
    ax1.plot(x_coordinates, y_coordinates, 'o',
                 color='w', marker='o', markersize=7, markeredgecolor='black')#, linewidth=0)

    # Seperate the stroke and control means and draw lines
    x_stroke = np.array([x_coordinates[0],x_coordinates[2],x_coordinates[4]])
    y_stroke = np.array([y_coordinates[0],y_coordinates[2],y_coordinates[4]])
    x_control = np.array([x_coordinates[1],x_coordinates[3],x_coordinates[5]])
    y_control = np.array([y_coordinates[1],y_coordinates[3],y_coordinates[5]])
    xp = np.linspace(0, 6.5, 10)

    slope, intercept, r_value, p_value, std_err = sp.stats.linregress(x_stroke,y_stroke)
    yp = xp*slope + intercept
    ax1.plot(xp,yp,color='#000000')
    print('stroke slope high freq: ',slope)

    slope, intercept, r_value, p_value, std_err = sp.stats.linregress(x_control,y_control)
    yp = xp*slope + intercept
    ax1.plot(xp,yp,color='#000000')
    print('control slope high freq: ',slope)

    fig.savefig('high_freq_content_boxplot.png')

    # frequency ratio box plot
    num_per_support = 3*5
    data_boxplot_low /= num_per_support

    figure_size = (8,3.55) # sets the size of the figure in inches

    fig_low, ax2 = plt.subplots(figsize=figure_size)
    print(data_boxplot_low)
    labels = ['Stroke-0%','Controls-0%','Stroke-20%','Controls-20%','Stroke-50%','Controls-50%']
    medianprops = dict(linewidth=2.5, color='black')
    bp = ax2.boxplot(data_boxplot_low, notch=0, medianprops=medianprops, labels=labels)#, patch_artist=True,boxprops=dict(facecolor=color_combine, color=c))
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='+')
    # Add a horizontal grid to the plot, but make it very light in color
    # so we can use it for reading data values but not be distracting
    ax2.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                   alpha=0.5)

    # Hide these grid behind plot objects
    ax2.set_axisbelow(True)
    ax2.set_title('Low-Frequency Content For Each Loading Level')#, fontsize=10, fontweight='bold')
    ax2.set_xlabel('Loading Level (% Max)')
    ax2.set_ylabel('Low-Frequency Content')
    # for label in (ax2.get_xticklabels() + ax2.get_yticklabels()):
    #     label.set_fontsize(8)

    box_colors = ['#601A4A', '#601A4A','#EE442F','#EE442F', '#63ACBE', '#63ACBE']
    box_alpha = [None,0.3,None,0.3,None,0.3]
    medians = np.empty(6)
    for i in range(6):
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        box_coords = np.column_stack([boxX, boxY])
        ax2.add_patch(Polygon(box_coords, facecolor=box_colors[i], alpha=box_alpha[i]))

    x_coordinates = np.array([1,2,3,4,5,6])
    y_coordinates = np.mean(data_boxplot_low,axis=0)
    # print(np.mean(data_boxplot_low,axis=0))
    ax2.plot(x_coordinates, y_coordinates, 'o',
                 color='w', marker='o', markersize=7, markeredgecolor='black')#, linewidth=0)

    # Seperate the stroke and control means and draw lines
    x_stroke = np.array([x_coordinates[0],x_coordinates[2],x_coordinates[4]])
    y_stroke = np.array([y_coordinates[0],y_coordinates[2],y_coordinates[4]])
    x_control = np.array([x_coordinates[1],x_coordinates[3],x_coordinates[5]])
    y_control = np.array([y_coordinates[1],y_coordinates[3],y_coordinates[5]])
    xp = np.linspace(0, 6.5, 10)

    slope, intercept, r_value, p_value, std_err = sp.stats.linregress(x_stroke,y_stroke)
    yp = xp*slope + intercept
    ax2.plot(xp,yp,color='#000000')
    print('slope low freq: ',slope)

    slope, intercept, r_value, p_value, std_err = sp.stats.linregress(x_control,y_control)
    yp = xp*slope + intercept
    ax2.plot(xp,yp,color='#000000')
    print('slope low freq: ',slope)

    fig_low.savefig('low_freq_content_boxplot.png')

plt.show()
