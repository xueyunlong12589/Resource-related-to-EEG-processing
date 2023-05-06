import os
import numpy as np
import math
import scipy.io as sio
from scipy.fftpack import fft
from scipy.signal.windows import hann,hanning


def DE_PSD(data,stft_para):
    '''
    compute DE and PSD
    --------
    input:  data [n*m]          n electrodes, m time points:62 5
            stft_para.stftn     frequency domain sampling rate:256
            stft_para.fStart    start frequency of each frequency band:[1,4,8,14,31]
            stft_para.fEnd      end frequency of each frequency band:[4,8,14,31,50]
            stft_para.window    window length of each sample point(seconds):1
            stft_para.overlap   overlap between adjacent samples(seconds):0
            stft_para.fs        original frequency:1000
    output:PSD,DE,PSD_MOV,DE_MOV [l*n*k]       l windows, n electrodes,  k frequency bands
    '''
    #initialize the parameters
    STFTN=stft_para['stftn']
    fStart=stft_para['fStart']
    fEnd=stft_para['fEnd']
    fs=stft_para['fs']
    window=stft_para['window']
    overlap=stft_para['overlap']

    fStartNum=np.zeros([len(fStart)],dtype=int)
    fEndNum=np.zeros([len(fEnd)],dtype=int)
    for i in range(0,len(fStart)):
        fStartNum[i]=int(fStart[i]/fs*STFTN)
        fEndNum[i]=int(fEnd[i]/fs*STFTN)

    #print(fStartNum[0],fEndNum[0])
    n=data.shape[0]
    m=data.shape[1]
    l=int((m-overlap)/(window-overlap))

    #print(m,n,l)
    PSD = np.zeros([n,len(fStart)])
    DE = np.zeros([n,len(fStart)])
    #Hanning window
    WindowPoints=fs*window
    #Hwindow=hanning(Hlength)
    Hwindow=hann(WindowPoints)
    
    for i in range(0,l):
        pstart=np.floor((i*(window-overlap))*fs)
        pend=np.floor(pstart+window*fs)
        dataNow=data[:,pstart:pend]
        for j in range(0,n):
            temp=dataNow[j]
            Hdata=temp*Hwindow
            FFTdata=fft(Hdata,STFTN)
            magFFTdata=abs(FFTdata[0:int(STFTN/2)])
            for p in range(0,len(fStart)):
                E = 0
                #E_log = 0
                for p0 in range(fStartNum[p]-1,fEndNum[p]):
                    E=E+magFFTdata[p0]*magFFTdata[p0] 
                #    E_log = E_log + log2(magFFTdata(p0)*magFFTdata(p0)+1)
                E = E/(fEndNum[p]-fStartNum[p]+1)
                PSD[i][j][p] = E
                DE[i][j][p] = math.log(100*E,2)
                #de(j,i,p)=log2((1+E)^4)
    
    PSD_MOV=smoothdata(PSD,5)
    DE_MOV=smoothdata(DE,5)
    
    return PSD,DE,PSD_MOV,DE_MOV

def smoothdata(data,window):
    #moving average
    out=np.zeros(data.shape)
    left=math.ceil((window-1)/2)
    right=math.floor((window-1)/2)
    n=data.shape[0]
    for i in range(n):
        pleft=max(0,i-left)
        pright=min(n-1,i+right)
        pdata=data[pleft:pright+1]
        out[i]=np.mean(pdata,0)
    return out


