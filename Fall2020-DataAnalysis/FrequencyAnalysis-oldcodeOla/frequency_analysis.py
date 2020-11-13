
import numpy as np
from numpy import fft, genfromtxt
from scipy import fft, arange
import matplotlib.pyplot as plt
import csv

# Edit these variables before running
DIR = "Z:" #set directory where data is mounted Ola- "/media/ola/Elements/R01prelim" Milli -"Z:"
number_of_control_subjects = 6
number_of_stroke_subjects = 6
DT = 0.01
Fs = 1/DT

def calculate_amplitude(w,y):
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
    T = n/Fs
    frq = k/T # two sides frequency range
    Y = abs(fft(y)) # fft computing (note that you might otherwise normalize by (2/n)
            # abs takes the magnitude of a complex number
    # Choose the one sided frequencies of interest
    frq = frq[range(n//2)] # one side frequency range
    Y = Y[range(n//2)]

    # Place values into prespecified frequency bins
    w_len = len(w)
    dw = w[1]-w[0]
    w_index = 0
    A = np.zeros(w_len)
    num = np.zeros(w_len)
    Y_len = len(Y)
    for i in range(0,Y_len):
        if frq[i]>(w[w_index]+(dw/2)):
            w_index += 1
        # print('i', i, 'frq[i]', frq[i], 'w_index', w_index, '(w[w_index]+(dw/2))', (w[w_index]+(dw/2)) )
        # print('w_index', w_index, 'w_len', w_len)
        if (w_index>=w_len):
            break
        A[w_index] += Y[i]
        num[w_index] += 1
        # passed = i/(Y_len-1)
        # print(passed)
    # print('here3')

    for i in range(w_len):
        A[i] = A[i]/num[i]

    E_signal = np.sum(np.square(A[1:w_len]))*(w[1]-w[0])
    # Calculate new signal normalized to have an energy of 1
    A *= np.sqrt(1/E_signal)
    E_signal = np.sum(np.square(A[1:w_len]))*(w[1]-w[0])
    print('Signal energy: ',E_signal) # should equal to 1

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
nyquist_freq = 15 # no valuable information was found after 15Hz
freq_step = 0.1
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
colors = ['r','b','g']
group_name = ['stroke','controls']

for group in range(len(group_name)):
    num_0percent = 0
    num_20percent = 0
    num_50percent = 0
    for subject_num in range(1,7):
        if (group_name[group] == 'stroke'):
            subname = "S0" + str(subject_num)
        else:
            subject_num += 10
            subname = "S" + str(subject_num)
        subfile = "/OutputData_S" + str(subject_num)
        plt.figure(subject_num)
        plt.title('Subject '+str(subject_num)+' Frequency Spectrum')
        for j in range(0,3): #support level
            num_supportlevel = 0
            for k in range(1,46): #trials
                for i in range(1,6): #Task number
                    try:
                        trialfile = DIR+subname+subfile+"_Task"+str(i)+"_SL"+str(j)+"_Trial"+ str(k)+".csv"
                        data = genfromtxt(trialfile,delimiter=',',dtype=float)
                        data = np.delete(data,0,0) # Deletes the first column of column names
                        print(trialfile)
                        delay = data[-1,0]-20.0 # Gets the final time to calculate the delay in data collection
                        delayInd = int(delay/DT)
                        score = data[-1,7]
                        trial_complete_index = int(data.shape[0]-1)
                        if score==20:
                            while data[trial_complete_index,7]==20:
                                trial_complete_index-=1
                        print('index: ',trial_complete_index, 'time: ', data[trial_complete_index,0]-delay)
                        # Combine the x and y direction forces
                        y = np.sqrt(np.square(data[delayInd:trial_complete_index,8])+np.square(data[delayInd:trial_complete_index,9]))
                        A = calculate_amplitude(w,y)

                        # For making individual plot
                        if (num_supportlevel==0):
                            A_support = A
                        else:
                            A_support += A
                        num_supportlevel += 1
                        # For combining stroke participants
                        if (j==0):
                            if (num_0percent==0):
                                A_0percent = A
                            else:
                                A_0percent += A
                            num_0percent += 1
                        if (j==1):
                            if (num_20percent==0):
                                A_20percent = A
                            else:
                                A_20percent += A
                            num_20percent += 1
                        if (j==2):
                            if (num_50percent==0):
                                A_50percent = A
                            else:
                                A_50percent += A
                            num_50percent += 1
                    except:
                        pass
            A_support /= num_supportlevel
            plt.plot(w[1:w_len],A_support[1:w_len],color=colors[j],label=supportlevels[j]) # plotting the spectrum (start at )
        plt.legend(loc="upper right")
        plt.xlabel('Freq (Hz)')
        plt.ylabel('Amplitude')
        plt.savefig('sub'+str(subject_num)+'_freq.png')
    # Group figure
    plt.figure(group+20)
    plt.title(group_name[group] + ' frequency spectrum')
    A_0percent /= num_0percent
    A_20percent /= num_20percent
    A_50percent /= num_50percent
    plt.plot(w[1:w_len],A_0percent[1:w_len],color=colors[0],label=supportlevels[0])
    plt.plot(w[1:w_len],A_20percent[1:w_len],color=colors[1],label=supportlevels[1])
    plt.plot(w[1:w_len],A_50percent[1:w_len],color=colors[2],label=supportlevels[2])
    plt.xlabel('Freq (Hz)')
    plt.ylabel('Amplitude')
    plt.legend(loc="upper right")
    plt.savefig(group_name[group]+'_freq.png')

    # Figure containing both stroke and controls
    plt.figure(30)
    plt.title('Frequency spectrum for stroke and controls')
    if (group_name[group] == 'stroke'):
        plt.plot(w[1:w_len],A_0percent[1:w_len],color=colors[0],label=('stroke '+supportlevels[0]))
        plt.plot(w[1:w_len],A_20percent[1:w_len],color=colors[1],label=('stroke '+supportlevels[1]))
        plt.plot(w[1:w_len],A_50percent[1:w_len],color=colors[2],label=('stroke '+supportlevels[2]))
    else:
        plt.plot(w[1:w_len],A_0percent[1:w_len],linestyle='--',color=colors[0],label=('control '+supportlevels[0]))
        plt.plot(w[1:w_len],A_20percent[1:w_len],linestyle='--',color=colors[1],label=('control '+supportlevels[1]))
        plt.plot(w[1:w_len],A_50percent[1:w_len],linestyle='--',color=colors[2],label=('control '+supportlevels[2]))

    plt.xlabel('Freq (Hz)')
    plt.ylabel('Amplitude')
    plt.legend(loc="upper right")
    plt.savefig('both_freq.png')


plt.show()
