import numpy as np
from numpy import fft
from numpy.fft import fft

def normalize_spectrum(A,dw):
    """
    Normalizes a frequency spectrum so that total energy =1
    Inputs:
        A - frequency amplitudes equally spaced
        dw - change in frequency between every two amplitudes

    Outputs:
        A - normalized frequency spectrum
    """
    # E_signal = np.sum(np.square(A))*dw # calculates energy
    E_signal = np.sum(np.square(A[1:]))*dw # calculates energy
    A *= np.sqrt(1/E_signal) # normalizes
    return A

def calculate_amplitude(w,y,Fs,type='bin'):
    """
    Calculates the discrete time fourier transform and frequency amplitudes for the signal,
    places the values in frequency "bins" specified by vector w or interpolates,
    and normalizes so the energy of the spectrum is 1
    Inputs:
        w - frequency values for "bins", note that it must start at 0
        y - discrete values of function
        Fs - the sampling frequency of the data
        type - either 'bin' for placing amplitudes into frequency ranges (necessary
                when comparing participants)
               'interp' for interpolates to get the amplitude for specific freuqencies
               or 'None'
    Outputs:
        A - vector of signal amplitude at each frequency in w (normalized)
        frq - the frequencies corresponding to the amplitudes (if binned,
                frq will be the same as w)
    """

    # Use FFT to calculate the frequency spectrum for the given data length
    n = len(y) # length of the signal
    Y = abs(np.fft.fft(y, n=n, norm='ortho')) #Normalize by signal length
    # Y = abs(np.fft.fft(y, n=n)) #No normalization
    frq = np.fft.fftfreq(n, d=1/Fs) # Determine the corresponding frequency values for the fft output
    # Choose the one sided frequencies of interest
    frq = frq[range(n//2)] # one side frequency range
    Y = Y[range(n//2)]
    # print(n,len(Y))
    dw = w[1]-w[0]
    # print('frq:', frq)
    # print('Y:', Y)
    # Place values into prespecified frequency bins, each with range (w[w_index]-(dw/2),w[w_index]+(dw/2)]
    if type=='bin':
        w_len = len(w)
        # dw = w[1]-w[0]
        # w_index = 0
        A = np.zeros(w_len)
        num = np.zeros(w_len)
        Y_len = len(Y)
        for i in range(0,Y_len): # iterate through the frequency outputs of fft
            for j in range(0,w_len):
                if frq[i]>=(w[j]-(dw/2)):
                    if frq[i]<(w[j]+(dw/2)):
                        A[j] += Y[i]
                        num[j] += 1

        # Take the average of each "bin"
        # A = np.divide(A,num)
        for i in range(w_len):
            if (num[i]>0):
                A[i] = A[i]/num[i]
            # else:
            #     A[i] = np.interp(w[i], frq, Y)
        frq = w

    elif type=='interp':
        A = np.interp(w, frq, Y)
        frq = w

    else:
        A = Y
        dw = frq[1]-frq[0]

    A = normalize_spectrum(A,dw)

    return A,frq
