function OutputSignal = denoise_harmonic_lf_auto(InputSignal)
%this code is designed to estimate and remove narrow bandwidth signals from
%data.  it is definitely good for LF data.  the user has to specify a
%vector of frequencies for removal.
%this one works by successively applying thresholding and narrow band
%filtering to estimate the signal in a narrow bandwidth, and then it
%subtracts the noise estimate

%written and edited by SAC July 2013

%this one is a good first one to use and only works above 9.1 kHz

%NOTES: working pretty good for LF.  need to handle lower frequencies, need to
%account for nearby in frequency noise level if wanted, and maybe consider
%a gaussian shape to the windowing instead of a boxcar to smooth out the
%edges of the removed noise

%make sure input signal is in InputSignal
%output signal will be in OutputSignal

%this one tries to go all the way to DC
%also decreases bandwidth down at lower frequencies
%also chooses noise level based on local background

%always subtract mean
InputSignal=InputSignal-mean(InputSignal);

%parameters to set
Fsamp=1e6; %sampling frequency, Hz

BWInit=9e3; %initial full bandwidth
BWHalvings(1)=6;
STDLevel=2.0; %thresholding level for saturating time domain signal
%for lower frequencies
BWHalvings(2)=10;


%inform the removal BW
FinalBW=BWInit/2^BWHalvings(1);
%disp(FinalBW);


%frequency of fft
freqvec=(1:length(InputSignal)-1)*Fsamp/length(InputSignal);
df=Fsamp/length(InputSignal);

OutputSignal=InputSignal;

KeepGoing=1;
FreqCount=0;
MaxFreqs=400;
MinFreq=9.1e3; %Hz

while(KeepGoing),
    FreqCount=FreqCount+1;
    
    %Compute current spectrum
    InputSpec=abs(fft(OutputSignal));
    
    %     %this section computes a threshold based on the entire signal
    %     [MaxVal,MaxLoc]=max(InputSpec(round(MinFreq/df):(length(InputSignal)/2)));
    %     Thresh2=7*std(InputSpec(round(MinFreq/df):(length(InputSignal)/2)));
    %     MaxLoc=MaxLoc+round(MinFreq/df)-1;
    
    %better might be
    %compute smoothed spectrum
    SmoothSpec=conv(InputSpec,ones(1,1001)/1001,'same');
    NormSpec=InputSpec./SmoothSpec;
    [MaxVal,MaxLoc]=max(NormSpec(round(MinFreq/df):(length(InputSignal)/2)));
    MaxLoc=MaxLoc+round(MinFreq/df)-1;
    Thresh2=5; %this is a normalized value that works well
    
    if MaxVal<Thresh2, %then no more noise frequencies
        KeepGoing=0;%disp('Ended normally.');
    end;
    if FreqCount>MaxFreqs,%then exceeded maximum
        KeepGoing=0;%disp('Reached MaxFreqs.');
    end;
    
    if KeepGoing,
        
        Fnoise=freqvec(MaxLoc);
        
        WorkingVec=OutputSignal;
        WorkingBW=BWInit;
        
        %find sample that is the closest to the center frequency
        [tmp,centersamp]=min(abs(freqvec-Fnoise));
        
        %choose bandwidth based on frequency
        if freqvec(centersamp)>0,
            SelectedHalvings=BWHalvings(1);
        else
            SelectedHalvings=BWHalvings(2);
        end;
        
        for k=1:SelectedHalvings,
            
            %threshold
            thresh=STDLevel*std(WorkingVec);
            WorkingVec(find(WorkingVec>thresh))=thresh;
            WorkingVec(find(WorkingVec<-thresh))=-thresh;
            
            %filter
            dataspec=fft(WorkingVec);
            dataspec2=zeros(size(dataspec));
            HalfWidthInSamples=round(WorkingBW/2/df);
            
            % win1=hamming(centersamp+HalfWidthInSamples).';
            % win1b=hamming(centersamp+HalfWidthInSamples-1).';
            % win2=hamming(length([-HalfWidthInSamples:HalfWidthInSamples])).';
            
            win1=rectwin(centersamp+HalfWidthInSamples).';
            win1b=rectwin(centersamp+HalfWidthInSamples-1).';
            win2=rectwin(length([-HalfWidthInSamples:HalfWidthInSamples])).';
            
            if centersamp-HalfWidthInSamples<=1,
                dataspec2(1:(centersamp+HalfWidthInSamples))=win1.*dataspec(1:(centersamp+HalfWidthInSamples));
                dataspec2(length(InputSignal)+2-centersamp+[-HalfWidthInSamples:(centersamp-2)])=win1b.*dataspec(length(InputSignal)+2-centersamp+[-HalfWidthInSamples:(centersamp-2)]);
            else
                dataspec2(centersamp+[-HalfWidthInSamples:HalfWidthInSamples])=win2.*dataspec(centersamp+[-HalfWidthInSamples:HalfWidthInSamples]);
                dataspec2(length(InputSignal)-centersamp+2+[-HalfWidthInSamples:HalfWidthInSamples])=win2.*dataspec(length(InputSignal)-centersamp+2+[-HalfWidthInSamples:HalfWidthInSamples]);
            end;
            
            WorkingVec=real(ifft(dataspec2));
            WorkingBW=round(WorkingBW/2);
        end;
        
        %plot(WorkingVec);pause;
        
        OutputSignal=OutputSignal-WorkingVec;
    end;
    
end;
% disp([num2str(FreqCount) ' frequencies cleaned.']);

%maybe do one last LPF if needed
%[b,a]=cheby1(4,0.1,2*300/1e3);
%OutputSignal=filtfilt(b,a,OutputSignal);



