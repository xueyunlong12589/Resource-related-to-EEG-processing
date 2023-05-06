function [ PSD, DE, PSD_MOV, DE_MOV ] = STFT(data,stft_para)
%input: data [n*m]          n electrodes, m time points:62 5
%       stft_para.stftn     frequency domain sampling rate:256
%       stft_para.fStart    start frequency of each frequency band:[1,4,8,14,31]
%       stft_para.fEnd      end frequency of each frequency band:[4,8,14,31,50]
%       stft_para.window    window length of each sample point(seconds):1
%       stft_para.overlap   overlap between adjacent samples(seconds):0
%       stft_para.f         original frequency:1000
%output:PSD,DE,PSD_MOV,DE_MOV [l*n*k]       l windows, n electrodes,  k frequency bands
   
    %initialize the parameters
    STFTN=stft_para.stftn;
    fStart=stft_para.fStart;
    fEnd=stft_para.fEnd;
    fs=stft_para.fs;
    window=stft_para.window;
    overlap=stft_para.overlap;
    
    fStartNum=zeros(1,length(fStart));
    fEndNum=zeros(1,length(fEnd));
    for i=1:length(fStart)
        fStartNum(1,i)=fix(fStart(i)/fs*STFTN);
        fEndNum(1,i)=fix(fEnd(i)/fs*STFTN);
    end
    
    [n,m]=size(data);
    l=fix((m-overlap)/(window-overlap));
    PSD=zeros(l,n,length(fStart));
    DE = zeros(l,n,length(fStart));
    %Hanning window
    WindowPoints=fs*window;
    Hwindow=hann(WindowPoints);
    
    for i=1:l
        pstart=floor(((i-1)*(window-overlap))*fs+1);
        pend=floor(pstart+window*fs-1);
        dataNow=data(:,pstart:pend);
        for j=1:n
            temp=dataNow(j,:);
            Hdata=temp.*Hwindow';
            FFTdata=fft(Hdata,STFTN);
            magFFTdata=abs(FFTdata(1:1:STFTN/2));
            for p=1:length(fStart)
                E=0;
                %E_log = 0;
                for p0=fStartNum(p):fEndNum(p)
                    E=E+magFFTdata(p0)*magFFTdata(p0);
                %    E_log = E_log + log2(magFFTdata(p0)*magFFTdata(p0)+1);
                end
                E=E/(fEndNum(p)-fStartNum(p)+1);
                PSD(i,j,p)= E;
                DE(i,j,p) = log2(100*E);
                %de(j,i,p)=log2((1+E)^4);
            end
        end
    end
    PSD_MOV=smoothdata(PSD,'movmean' ,5);
    DE_MOV=smoothdata(DE,'movmean' ,5);
end


