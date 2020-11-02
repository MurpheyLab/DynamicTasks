import numpy as np
from numpy import fft

# This version does not seperate into 'bins'
def calculate_amplitude(w,y,Fs,bin):
    """
    Calculates the discrete time fourier transform and frequency amplitudes for the signal,
    places the values in frequency "bins" specified by vector w,
    and normalizes so the energy of the spectrum is 1
    Inputs:
        w - frequency values for "bins", note that it must start at 0
        y - discrete values of function
        Fs - the sampling frequency of the data
        bin - boolean: True - bin amplitudes into frequency ranges (necessary
                when comparing participants), False - do not bin
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
    dw = frq[1]-frq[0]
    # print('frq:', frq)
    # print('Y:', Y)
    # Place values into prespecified frequency bins, each with range (w[w_index]-(dw/2),w[w_index]+(dw/2)]
    if bin:
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
        frq = w
    else:
        A = Y

    # Calculate the energy of the amplitude signal
    E_signal = np.sum(np.square(A))*dw
    # print('Signal energy: ',E_signal)

    # Calculate new signal normalized to have an energy of 1
    A *= np.sqrt(1/E_signal)
    E_signal = np.sum(np.square(A))*dw
    # print('Signal energy: ',E_signal) # should equal to 1

    return A,frq
