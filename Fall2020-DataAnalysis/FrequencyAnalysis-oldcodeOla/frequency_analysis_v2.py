
import numpy as np
import scipy as sp
from numpy import fft, genfromtxt
from scipy import fft, arange, signal
import matplotlib.pyplot as plt
import csv

# Edit these variables before running
DIR = "/media/ola/Data/Research/DowntownData/Fall2020LabData/"
number_of_control_subjects = 6
number_of_stroke_subjects = 6
DT = 0.05
Fs = 1/DT
freq_timestep = DT

def frequency_analysis(y):
    """
    Calculates the discrete time fourier transform and frequency amplitudes for the signal,
    places the values in frequency bins specified by vector w,
    and normalizes so the energy of the spectrum is 1
    Inputs:
        w - frequency values
        y - discrete values of function
    Outputs:
        A - vector of signal amplitude at each frequency in w (normalized)
    """

    # Use FFT to calculate the frequency spectrum for the given data length
    n = len(y) # length of the signal
    k = arange(n)

    signal = y
    fourier = abs(np.fft.fft(signal, norm='ortho'))
    # print("fourier:", len(fourier))
    # print(fourier)
    freq = np.fft.fftfreq(len(fourier), d=freq_timestep)

    E_signal = np.sum(np.square(fourier[:]))*(w[1]-w[0])
    # Calculate new signal normalized to have an energy of 1
    fourier *= np.sqrt(1/E_signal)
    E_signal = np.sum(np.square(fourier[:]))*(w[1]-w[0])
    #print('Signal energy: ',E_signal) # should equal to 1

    # find local maxima
    peaks = (sp.signal.argrelmax(fourier[:len(freq)/2], order=12, mode='wrap')[0]).tolist()
    # print(peaks)
    peak_indicies = peaks
    peak_amplitude = fourier[peaks]

    # find tallest local maxima
    num_of_max = 3
    tallest_peaks = np.zeros((num_of_max,2))
    for i in range(num_of_max):
        arg = np.argmax(peak_amplitude)
        tallest_peaks[i,0] = peak_indicies[arg]
        tallest_peaks[i,1] = peak_amplitude[arg]
        peak_indicies = np.delete(peak_indicies,arg)
        peak_amplitude = np.delete(peak_amplitude,arg)
    # print('tallest peaks: ', tallest_peaks)

    # plt.figure(1)
    # plt.plot(freq[:len(freq)/2],fourier[:len(freq)/2])
    # #peaks2 = [peaks[i]+len(freq)/2 for i in range(len(freq)/2)]
    # plt.plot(freq[peaks], fourier[peaks], "x")
    # ind = [int(tallest_peaks[i,0]) for i in range(len(tallest_peaks))]
    # plt.plot(freq[ind], tallest_peaks[:,1], "x")
    # plt.show()
    # print(fourier)
    return [fourier, freq, tallest_peaks]


# Create vector of frequencies of interest
nyquist_freq = 10 # no valuable information was found after 15Hz
freq_step = 0.1
w = np.zeros(int(np.floor(nyquist_freq/freq_step)))
frq = 0 # freq_step
for i in range(len(w)):
    w[i] = frq
    frq = frq + freq_step
#w[-1] = frq

print('w: ', w)
w_len = len(w)
# Label factors as strings for ezANOVA analysis
tasks = ['Task1','Task2','Task3','Task4','Task5']
subjects = ['Subject2','Subject3','Subject4','Subject5','Subject6','Subject7']
colors = ['r','b','g','k']


total = 0
num_freq_nf = [0,0,0,0]
A_freq_nf = [[],[],[],[]]
num_freq_wf = [0,0,0,0]
A_freq_wf = [[],[],[],[]]
for subject_num in range(2,8):
    print(subject_num)
    subname = "S0" + str(subject_num)
    subfile = "/OutputData_S" + str(subject_num)
    for j in range(0,4): #frequency
        for k in range(1,49): #trials
            for i in range(0,2): #forces
                try:
                    trialfile = DIR+subname+subfile+"_Trial"+str(k)+"_Freq"+str(j)+"_SL0_F"+str(i)+".csv"
                    data = genfromtxt(trialfile,delimiter=',',dtype=float)
                    total += 1
                    data = np.delete(data,0,0) # Deletes the first column of column names
                    # print(trialfile)
                    score = data[-1,7]
                    # Combine the x and y direction forces
                    # y = np.sqrt(np.square(data[:590,8])+np.square(data[:590,9]))
                    #y = np.square(data[:590,8]) # just x
                    y = np.square(data[:590,9]) # just y
                    [A, bins, tallpeaks] = frequency_analysis(y)
                    # For combining responses at the same freq
                    if (i==0):
                        if num_freq_nf[j]==0:
                            A_freq_nf[j] = A
                        else:
                            A_freq_nf[j] += A
                        num_freq_nf[j] += 1
                    if (i==1):
                        if num_freq_wf[j]==0:
                            A_freq_wf[j] = A
                        else:
                            A_freq_wf[j] += A
                        num_freq_wf[j] += 1

                    ###### write data to a CSV file to do statistics on
                    # TO DO --------------
                except:
                    pass
    #print(total, 48*6)

    ##### individual plots
    # plt.figure(subject_num)
    # for j in range(0,4):
    #     plt.plot(bins[:len(bins)/2],A_freq_nf[j][:len(bins)/2],color=colors[j],linestyle='dashed',label="Freq"+str(j)+" w/o forces")
    #     plt.plot(bins[:len(bins)/2],A_freq_wf[j][:len(bins)/2],color=colors[j],label="Freq"+str(j)+" w/ forces") # plotting the spectrum (start at )
    # plt.legend(loc="upper right")
    # plt.xlabel('Freq (Hz)')
    # plt.ylabel('Amplitude')
    # plt.ylim(0,50)
    # plt.xlim(0,10)
    # plt.savefig(DIR+subname+'/sub'+str(subject_num)+'_freq.png')


def shrink(data, rows, cols):
    return data.reshape(rows, data.shape[0]/rows, cols, data.shape[1]/cols).sum(axis=1).sum(axis=2)

for j in range(0,4):
    plt.figure(10+j)
    plt.plot(bins[:len(bins)/2],A_freq_nf[j][:len(bins)/2]/num_freq_nf[j],color=colors[j],linestyle='dashed',label="Freq"+str(j)+" w/o forces")
    plt.plot(bins[:len(bins)/2],A_freq_wf[j][:len(bins)/2]/num_freq_wf[j],color=colors[j],label="Freq"+str(j)+" w/ forces")
    #plt.plot(bins[:len(bins)/2],A_freq_wf[j][:len(bins)/2]-A_freq_nf[j][:len(bins)/2],color=colors[j],label="Freq"+str(j)+" w/ forces")
    plt.legend(loc="upper right")
    plt.xlabel('Freq (Hz)')
    plt.ylabel('Amplitude')
    #plt.ylim(0,10)
    #plt.xlim(0,10)
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig(DIR+'aggregate_freq'+str(j)+'.png')

    ############ Plotting ################
    #
    # ###### Group figure
    # plt.figure(group+20)
    # plt.title(group_name[group] + ' frequency spectrum')
    # A_0percent /= num_0percent
    # A_20percent /= num_20percent
    # A_50percent /= num_50percent
    # plt.plot(bins[:len(bins)/2],A_0percent[:len(bins)/2],color=colors[0],label=supportlevels[0])
    # plt.plot(bins[:len(bins)/2],A_20percent[:len(bins)/2],color=colors[1],label=supportlevels[1])
    # plt.plot(bins[:len(bins)/2],A_50percent[:len(bins)/2],color=colors[2],label=supportlevels[2])
    # # plt.plot(w[1:w_len],A_0percent[1:w_len],color=colors[0],label=supportlevels[0])
    # # plt.plot(w[1:w_len],A_20percent[1:w_len],color=colors[1],label=supportlevels[1])
    # # plt.plot(w[1:w_len],A_50percent[1:w_len],color=colors[2],label=supportlevels[2])
    # plt.xlabel('Freq (Hz)')
    # plt.ylabel('Amplitude')
    # # plt.xlim(0,2)
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.legend(loc="upper right")
    # plt.savefig(group_name[group]+'_freq.png')
    #
    # ###### Figure containing both stroke and controls
    # plt.figure(30)
    # plt.title('Frequency spectrum for stroke and controls')
    # if (group_name[group] == 'stroke'):
    #     # plt.plot(w[1:w_len],A_0percent[1:w_len],color=colors[0],label=('stroke '+supportlevels[0]))
    #     # plt.plot(w[1:w_len],A_20percent[1:w_len],color=colors[1],label=('stroke '+supportlevels[1]))
    #     # plt.plot(w[1:w_len],A_50percent[1:w_len],color=colors[2],label=('stroke '+supportlevels[2]))
    #     plt.plot(bins[:len(bins)/2],A_0percent[:len(bins)/2],color=colors[0],label=('stroke '+supportlevels[0]))
    #     plt.plot(bins[:len(bins)/2],A_20percent[:len(bins)/2],color=colors[1],label=('stroke '+supportlevels[1]))
    #     plt.plot(bins[:len(bins)/2],A_50percent[:len(bins)/2],color=colors[2],label=('stroke '+supportlevels[2]))
    # else:
    #     plt.plot(bins[:len(bins)/2],A_0percent[:len(bins)/2],linestyle='--',color=colors[0],label=('control '+supportlevels[0]))
    #     plt.plot(bins[:len(bins)/2],A_20percent[:len(bins)/2],linestyle='--',color=colors[1],label=('control '+supportlevels[1]))
    #     plt.plot(bins[:len(bins)/2],A_50percent[:len(bins)/2],linestyle='--',color=colors[2],label=('control '+supportlevels[2]))
    #     # plt.plot(w[1:w_len],A_0percent[1:w_len],linestyle='--',color=colors[0],label=('control '+supportlevels[0]))
    #     # plt.plot(w[1:w_len],A_20percent[1:w_len],linestyle='--',color=colors[1],label=('control '+supportlevels[1]))
    #     # plt.plot(w[1:w_len],A_50percent[1:w_len],linestyle='--',color=colors[2],label=('control '+supportlevels[2]))
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.xlabel('Freq (Hz)')
    # plt.ylabel('Amplitude')
    # # plt.xlim(0,2)
    # plt.legend(loc="upper right")
    # plt.savefig('both_freq.png')


plt.show()
