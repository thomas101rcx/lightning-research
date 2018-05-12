%% Spectral analysis


% Choose which flash to analyze
[filename, pathname] = uigetfile({'*.mat';'*.*'},'File Selector');
sm1 = load([pathname,filename]);

% Choose the corresponding VHF data -> Convert fig to array;

[filename, pathname] = uigetfile({'*.fig'},'File Selector');
S2 = open([pathname,filename]);
h = gcf;
set(h, 'Visible', 'off');
axesObjs = get(h, 'Children');
dataObjs = get(axesObjs, 'Children');
time_temp = get(dataObjs{3}, 'XData');
VHF_time = [time_temp{2}]';
Azimuth_temp = get(dataObjs{3}, 'YData');
VHF_Azimuth = [Azimuth_temp{2}]';
Elevation_temp = get(dataObjs{4}, 'YData');
VHF_Elevation = [Elevation_temp{2}]';

S2 = [VHF_time VHF_Azimuth VHF_Elevation];

Flash_start_time = round(sm1.Duke_time(1),6);
First_Pulse_time = round(S2(1,1), 6);
Last_Pulse_time = round(S2(length(S2(:,1)),1),6);
Max_Interval_VHF = Last_Pulse_time - First_Pulse_time;
First_Pulse_index = (First_Pulse_time - Flash_start_time)*10^(6);

Starttime = round(First_Pulse_time,4);
Interval = 0.02;
Endtime = Starttime + Interval;
globalstart = int32(round((Starttime - Flash_start_time)*10^(6),6));
globalend =  int32(globalstart + Interval*10^(6));

slidewindow = 100; %us
overlap = 25;  %us

% Transpose data

sig_Duke_BGE_unfiltered = sm1.Duke_BGE';
sig_Duke_BGN_unfiltered = sm1.Duke_BGN';

sig_Hudson_BGE_unfiltered = sm1.Hudson_BGE';
sig_Hudson_BGN_unfiltered = sm1.Hudson_BGN';

sig_PS2_BGE_unfiltered = sm1.PS2_BGE';
sig_PS2_BGN_unfiltered = sm1.PS2_BGN';

sig_PS3_BGE_unfiltered = sm1.PS3_BGE';
sig_PS3_BGN_unfiltered = sm1.PS3_BGN';

% Preflash standard deviation

Duke_BGE_preflashstd = single(std(sm1.Duke_BGE(1:1:int32(round((First_Pulse_time-Flash_start_time)*10^(6), 6)))));
Duke_BGN_preflashstd = single(std(sm1.Duke_BGN(1:1:int32(round((First_Pulse_time-Flash_start_time)*10^(6), 6)))));

Hudson_BGE_preflashstd = single(std(sm1.Hudson_BGE(1:1:int32(round((First_Pulse_time-Flash_start_time)*10^(6), 6)))));
Hudson_BGN_preflashstd = single(std(sm1.Hudson_BGN(1:1:int32(round((First_Pulse_time-Flash_start_time)*10^(6), 6)))));

PS2_BGE_preflashstd = single(std(sm1.PS2_BGE(1:1:int32(round((First_Pulse_time-Flash_start_time)*10^(6), 6)))));
PS2_BGN_preflashstd = single(std(sm1.PS2_BGN(1:1:int32(round((First_Pulse_time-Flash_start_time)*10^(6), 6)))));

PS3_BGE_preflashstd = single(std(sm1.PS3_BGE(1:1:int32(round((First_Pulse_time-Flash_start_time)*10^(6), 6)))));
PS3_BGN_preflashstd = single(std(sm1.PS3_BGN(1:1:int32(round((First_Pulse_time-Flash_start_time)*10^(6), 6)))));



% Interpolation 

sig_Duke_BGE = single(interp(sig_Duke_BGE_unfiltered(globalstart:1:globalend),10));
sig_Duke_BGN = single(interp(sig_Duke_BGN_unfiltered(globalstart:1:globalend),10));

sig_Hudson_BGE = single(interp(sig_Hudson_BGE_unfiltered(globalstart:1:globalend),10));
sig_Hudson_BGN = single(interp(sig_Hudson_BGN_unfiltered(globalstart:1:globalend),10));

sig_PS2_BGE = single(interp(sig_PS2_BGE_unfiltered(globalstart:1:globalend),10));
sig_PS2_BGN = single(interp(sig_PS2_BGN_unfiltered(globalstart:1:globalend),10));

sig_PS3_BGE = single(interp(sig_PS3_BGE_unfiltered(globalstart:1:globalend),10));
sig_PS3_BGN = single(interp(sig_PS3_BGN_unfiltered(globalstart:1:globalend),10));

PS3_Freq = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');

Totalpoints = int32((globalend - globalstart) / overlap) - 1;

% Spectral analysis

Flashname = filename(1:1:19);
directoryname = sprintf('%s_Flash_%dus_slidingwindow_%dus_overlap_%9.6f-%9.6f_spectral_analysis', Flashname,slidewindow, overlap, Starttime, Endtime);
mkdir(directoryname);
cd(directoryname);


parfor ( i = 1:1:Totalpoints-4, 4)
    
    
   [PS3_freq] =  spectral((i - 1) * overlap*10 + 1, (i - 1) * overlap*10 + slidewindow*10 ,sig_Duke_BGE ,sig_Duke_BGN,sig_Hudson_BGE,sig_Hudson_BGN, sig_PS2_BGE, sig_PS2_BGN, sig_PS3_BGE,sig_PS3_BGN , Duke_BGE_preflashstd, Duke_BGN_preflashstd, Hudson_BGE_preflashstd, Hudson_BGN_preflashstd, PS2_BGE_preflashstd, PS2_BGN_preflashstd, PS3_BGE_preflashstd, PS3_BGN_preflashstd, First_Pulse_time)
    
   PS3_Freq(i) = PS3_freq;
end

filename = sprintf('%dus_slidingwindow_%dus_overlap_%9.6f-%9.6f_All_window.mat', slidewindow, overlap, Starttime, Endtime);



function [PS3_freq] = spectral(startposition, Endposition, sig_Duke_BGE ,sig_Duke_BGN,sig_Hudson_BGE,sig_Hudson_BGN, sig_PS2_BGE, sig_PS2_BGN, sig_PS3_BGE,sig_PS3_BGN, Duke_BGE_preflashstd, Duke_BGN_preflashstd, Hudson_BGE_preflashstd, Hudson_BGN_preflashstd, PS2_BGE_preflashstd, PS2_BGN_preflashstd, PS3_BGE_preflashstd, PS3_BGN_preflashstd, First_Pulse_time)

%% Pre signal analysis ->
% 1. Filter out the low SNR signal using Max(Sig)/std(noise)
% 2. 
    
%     Determine which signal to use BGE or BGN based on the max(abs) value of the signal
%
%     Choose Duke as signal polarity reference
%
%     Signal 1 = BGE;
%     Signal 2 = BGN;

    SNR_Duke_BGE = max(abs(sig_Duke_BGE(startposition:1:Endposition)))/Duke_BGE_preflashstd;
    SNR_Duke_BGN = max(abs(sig_Duke_BGN(startposition:1:Endposition)))/Duke_BGN_preflashstd;
    
    SNR_Hudson_BGE = max(abs(sig_Hudson_BGE(startposition:1:Endposition)))/Hudson_BGE_preflashstd;
    SNR_Hudson_BGN = max(abs(sig_Hudson_BGN(startposition:1:Endposition)))/Hudson_BGN_preflashstd;
    
    SNR_PS2_BGE = max(abs(sig_PS2_BGE(startposition:1:Endposition)))/PS2_BGE_preflashstd;
    SNR_PS2_BGN = max(abs(sig_PS2_BGN(startposition:1:Endposition)))/PS2_BGN_preflashstd;

    SNR_PS3_BGE = max(abs(sig_PS3_BGE(startposition:1:Endposition)))/PS3_BGE_preflashstd;
    SNR_PS3_BGN = max(abs(sig_PS3_BGN(startposition:1:Endposition)))/PS3_BGN_preflashstd;
    
    Duke_signal_matrix_max = [SNR_Duke_BGE, SNR_Duke_BGN];
    Hudson_signal_matrix_max = [SNR_Hudson_BGE, SNR_Hudson_BGN];
    PS2_signal_matrix_max = [SNR_PS2_BGE, SNR_PS2_BGN];
    PS3_signal_matrix_max = [SNR_PS3_BGE, SNR_PS3_BGN];

    [~, Duke_signal] = max(Duke_signal_matrix_max);
    [~, Hudson_signal] = max(Hudson_signal_matrix_max);
    [~, PS2_signal] = max(PS2_signal_matrix_max);
    [~, PS3_signal] = max(PS3_signal_matrix_max);
    
    %% Select a higher SNR ratio to do Xcorr , either BGE or BGN 

    
   if (Duke_signal == 1)
        sig_Duke  = sig_Duke_BGE(startposition:1:Endposition);
        SNR_Duke = SNR_Duke_BGE;
        Max_Duke = max(abs(sig_Duke_BGE(startposition:1:Endposition)))/Duke_BGE_preflashstd;
    else
        sig_Duke  = sig_Duke_BGN(startposition:1:Endposition);
        SNR_Duke = SNR_Duke_BGN;
        Max_Duke = max(abs(sig_Duke_BGN(startposition:1:Endposition)))/Duke_BGN_preflashstd;
    end

    if (Hudson_signal == 1)
        sig_Hudson  = sig_Hudson_BGE(startposition:1:Endposition);
        SNR_Hudson = SNR_Hudson_BGE;
        Max_Hudson = max(abs(sig_Hudson_BGE(startposition:1:Endposition)));
    else
        sig_Hudson  = sig_Hudson_BGN(startposition:1:Endposition);
        SNR_Hudson = SNR_Hudson_BGN;
        Max_Hudson = max(abs(sig_Hudson_BGN(startposition:1:Endposition)));
    end

    if (PS2_signal == 1)
        sig_PS2  = sig_PS2_BGE(startposition:1:Endposition);
        SNR_PS2 = SNR_PS2_BGE;
        Max_PS2 = max(abs(sig_PS2_BGE(startposition:1:Endposition)))/PS2_BGE_preflashstd;
    else
        sig_PS2  = sig_PS2_BGN(startposition:1:Endposition);
        SNR_PS2 = SNR_PS2_BGN;
        Max_PS2 = max(abs(sig_PS2_BGN(startposition:1:Endposition)))/PS2_BGN_preflashstd;

    end

    if (PS3_signal == 1)
        sig_PS3  = sig_PS3_BGE(startposition:1:Endposition);
        SNR_PS3 = SNR_PS3_BGE;
        Max_PS3 = max(abs(sig_PS3_BGE(startposition:1:Endposition)))/PS3_BGE_preflashstd;

    else
        sig_PS3  = sig_PS3_BGN(startposition:1:Endposition);
        SNR_PS3 = SNR_PS3_BGN;
        Max_PS3 = max(abs(sig_PS3_BGN(startposition:1:Endposition)))/PS3_BGN_preflashstd;

    end
    [~, Duke_max_index] = max(sig_Duke);
    [~, Duke_min_index] = min(sig_Duke);
    [~, Hudson_max_index] = max(sig_Hudson);
    [~, Hudson_min_index] = min(sig_Hudson);
    [~, PS2_max_index] = max(sig_PS2);
    [~, PS2_min_index] = min(sig_PS2);
    [~, PS3_max_index] = max(sig_PS3);
    [~, PS3_min_index] = min(sig_PS3);
    
    % Find the max/min signal's index and return it 

    if Duke_max_index < Duke_min_index

            % The polarity of first pulse of Duke is positive

            if Hudson_max_index < Hudson_min_index  %Hudson is positive
                sig_Hudson = 1*sig_Hudson;
            else
                sig_Hudson = -1*sig_Hudson;
            end

            if PS2_max_index < PS2_min_index
                sig_PS2 = 1*sig_PS2;
            else
                sig_PS2 = -1*sig_PS2;
            end

            if PS3_max_index < PS3_min_index
                sig_PS3 = 1*sig_PS3;
            else
                sig_PS3 = -1*sig_PS3;
            end

    else

            % The polarity of first pulse of Duke is negative

            if Hudson_max_index > Hudson_min_index  %Hudson is negative
                sig_Hudson = 1*sig_Hudson;
            else
                sig_Hudson = -1*sig_Hudson;
            end

            if PS2_max_index > PS2_min_index
                sig_PS2 = 1*sig_PS2;
            else
                sig_PS2 = -1*sig_PS2;
            end

            if PS3_max_index > PS3_min_index
                sig_PS3 = 1*sig_PS3;
            else
                sig_PS3 = -1*sig_PS3;
            end

    end
    
    Max_std = [SNR_Duke,  SNR_Hudson,  SNR_PS2, SNR_PS3];

    
    f1 = figure(3);
    subplot(3,3,[1 2 3]);
    set(f1, 'visible', 'off');
    b = plot(1:1:length(sig_Duke),sig_Duke,'k-' , 1:1:length(sig_Duke), sig_Hudson, 'r', 1:1:length(sig_Duke), sig_PS2, 'b', 1:1:length(sig_Duke), sig_PS3, 'g');
    legend('Duke', 'Hudson', 'Ps2', 'Ps3');
    set(b(1), 'linewidth', 2);
    set(b(2), 'linewidth', 0.75);
    set(b(3), 'linewidth', 0.75);
    set(b(4), 'linewidth', 0.75);
    xlabel('time (0.1 us)');
    ylabel('BGE/BGN Magnitude');
    Abstimestart = First_Pulse_time + double(startposition-1)*0.0000001;
    Abstimeend = First_Pulse_time + double(Endposition-1)*0.0000001;
    title(sprintf('Original signal characteristics %9.6fus to %9.6fus',Abstimestart, Abstimeend), 'Fontsize', 9);
    values = sprintf('%5.2f ', Max_std);
    dim = [0.7 0.7 0.3 0.3];
    annotation('textbox',dim , 'String', values);
    
    Frametime = double(startposition)*0.0000001 + double(First_Pulse_time);
    
    fs = 10^(7);
    n = length(sig_PS3);
    fftPS3 = fft(sig_PS3);
    fftPS3_shift = fftshift(fftPS3);
    fshift = (-n/2:n/2-1)*(fs/n);
    powershift_PS3 = abs(fftPS3_shift).^2/n;
    maxfreq = fshift(powershift_PS3 == max(powershift_PS3));
    PS3_freq = max(maxfreq);
    

    fftDuke = fft(sig_Duke);
    fftDuke_shift = fftshift(fftDuke);
    powershift_Duke = abs(fftDuke_shift).^2/n;
    
    subplot(3, 3, [4 5 6]);
    c = plot(fshift, powershift_PS3, fshift, powershift_Duke);
    legend('PS3', 'Duke')
    
    
    d = designfilt('bandstopiir','FilterOrder',6, ...
               'HalfPowerFrequency1',310000,'HalfPowerFrequency2',330000, ...
               'DesignMethod','butter','SampleRate',fs);
    d1 = designfilt('bandstopiir','FilterOrder',6, ...
               'HalfPowerFrequency1',140000,'HalfPowerFrequency2',160000, ...
               'DesignMethod','butter','SampleRate',fs);
           
           
    filtered_PS3_Step1 = filtfilt(d, double(sig_PS3));
    filtered_PS3_final = filtfilt(d1, double(filtered_PS3_Step1));
    
    subplot(3, 3, [7 8 9]);
    e = plot(1:1:length(sig_Duke),sig_Duke,'k-' , 1:1:length(sig_Duke), sig_Hudson, 'r', 1:1:length(sig_Duke), sig_PS2, 'b', 1:1:length(sig_Duke), filtered_PS3_final, 'g');
    set(e(1), 'linewidth', 2);
    set(e(2), 'linewidth', 0.75);
    set(e(3), 'linewidth', 0.75);
    set(e(4), 'linewidth', 0.75);
    xlabel('time (0.1 us)');
    ylabel('BGE/BGN Magnitude');
    Abstimestart = First_Pulse_time + double(startposition-1)*0.0000001;
    Abstimeend = First_Pulse_time + double(Endposition-1)*0.0000001;
    title(sprintf('Original signal characteristics %9.6fus to %9.6fus',Abstimestart, Abstimeend), 'Fontsize', 9);
    values = sprintf('%5.2f ', Max_std);

    Frametime = double(startposition)*0.0000001 + double(First_Pulse_time);
    
    saveas(f1,sprintf('%f time.jpeg', Frametime), 'jpeg');
    delete(f1);
    
    
    
    clf
    




end
