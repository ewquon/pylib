#!/usr/bin/env python3
#
# Calculates wave spectrum for regular and irregular waves
# Based on Amy_PSD2.m from Fabian Wendt
#
# Written by Eliot Quon (eliot.quon@nrel.gov) -- 2017-12-12
#
import numpy as np
import matplotlib.pyplot as plt

def first_order_spectrum(time, x, smoothing=False, smoothtype='averaging', plotting=False, **kwargs):
    """Calculates the PSD or "first-order spectrum"

    Parameters
    ----------
    time : np.ndarray
        Time vector
    x : np.ndarray
        Data to be analyzed (vector)
    smoothing : boolean
        Flag for whether you want to perform smoothing. Note that you
        might need to play with parameters to get good smoothing.
    plotting : boolean
        Flag for whether you want to have results plotted
    smoothtype : string
        Options include 'averaging' (set parameters directly in code),
        'local' for localized smoothing (e.g., for irregular waves),
        'gaussian' for a Gaussian window
    **kwargs : dict
        Additional parameters, e.g. ni=2048 for averaging or nskip=5 
        for regular wave smoothing

    Outputs
    -------
    f : np.ndarray
        Vector of frequency values
    x_fft : np.ndarray
        One-sided FFT of x
    x_psd : np.ndarray
        One-sided PSD of x
    """

    """
    Calculate the PSD
    """
    N = len(time);
    dt = time[1] - time[0]
    time_span = time[-1] - time[0]
    Nfft = int(N/2)

    df = 1.0 / time_span;
    domega = df*(2*np.pi)
    f = np.arange(N,dtype=float)/(N-1) * (1.0/dt)  # Note Andy starts at df, rather than 0
    omega = 2*np.pi*f

    # Fourier tranform of signal
    x_fft = np.fft.fft(x,N);   

    # Calculate one-sided PSD (only use first half of resulting vector)
    x_psd = 2*np.abs( x_fft*np.conj(x_fft) )/(N*N*df);
    x_psd = x_psd[:Nfft+1]  # To get just one-sided

    """
    Smooth the PSD
    """
    if smoothing:
        if smoothtype == 'averaging':
            # first through averaging
            ni = kwargs['ni']
            nens = int(N / ni)
            x_psd2 = np.zeros((ni,))
            df2 = 1.0 / (time[ni]-time[0])
            f2 = np.arange(ni,dtype=float)/(ni-1) * (1.0/dt)

            for i in range(nens):
                x_fft2 = np.fft.fft.fft(x[ni*i:ni*(i+1)+1])
                x_psd2 = x_psd2 + 2*np.abs( x_fft2*np.conj(x_fft2) )/(ni*ni*df2)
            x_psd2 = x_psd2 / nens;
            x_psd2 = x_psd2[:ni/2+1]

        elif smoothtype == 'local':
            # then through MARIN's method
            nskip = kwargs['nskip']
            window = np.ones((nskip,)) / nskip
            x_fft2 = x_fft
            f2 = f
            #x_psd2 = convn(x_psd,window,'same');
            x_psd2 = np.convolve(x_psd,window,'same');

        elif smoothtype == 'gaussian':
            import scipy.signal
            nskip = kwargs['nskip']
            window = scipy.signal.gaussian(nskip)
            window = window / np.sum(window);
            x_fft2 = x_fft
            f2 = f
            #x_psd2 = conv(x_psd,window,'same');
            x_psd2 = np.convolve(x_psd,window,'same');
        
    """Plot PSDs"""
    if plotting:
        # 1-sided PSD
        fig,ax = plt.subplots(nrows=3)
        ax[0].plot(time,x)
        ax[0].set_xlabel('Time (sec)')
        ax[0].set_ylabel('Response');
        ax[1].plot(f[:Nfft+1],x_fft[:Nfft+1])
        ax[1].set_xlabel('Frequency (Hz)')
        ax[1].set_ylabel('Magnitude')
        ax[2].plot(f[:Nfft+1],x_psd[:Nfft+1])
        ax[2].set_xlabel('Frequency (Hz)')
        ax[2].set_ylabel('PSD')

        # 1-sided PSD on semilog plot
        fig,ax = plt.subplots()
        plt.semilogy(f[:Nfft+1],x_psd[:Nfft+1])
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('PSD (m^2/Hz)')
        plt.grid(True)
        plt.axis((0,1,10e-6,10e3))

        if smoothing:
            # plot smoothed PSD on top
            plt.semilogy(f2[:Nfft+1],x_psd2[:Nfft+1],'r')
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('PSD (m^2/Hz)')
            plt.grid(True)
            plt.axis((0,1,10e-6,10e3))
            plt.legend(['PSD','Averaged PSD'])

    # Output smoothed PSD if requested
    if smoothing:
        N = len(f2);
        x_psd = x_psd2
        f = f2[:Nfft+1]
        x_fft = x_fft2
    else:
        f = f[:Nfft+1]

    return f, x_fft, x_psd

#===============================================================================
if __name__ == '__main__':
    import sys
    fname = sys.argv[1]
    time,signal = np.loadtxt(fname,unpack=True)
    f0,fft0,psd0 = first_order_spectrum(time,signal,plotting=False,smoothing=False)
    f, fft, psd  = first_order_spectrum(time,signal,plotting=False,
                                        smoothing=True,smoothtype='local',nskip=2)

    plt.semilogy(f0,psd0,label='plain signal')
    plt.semilogy(f,psd,label='with smoothing')
    plt.xlim((-.05,.95))
    plt.ylim((5e-6,2e2))

    plt.show()
