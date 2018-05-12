%% Cross-corr algorithm for lightning detection
tic
clear; format short e

% mex linearto3d.c

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

Starttime = First_Pulse_time;
Interval = round(max(VHF_time)- min(VHF_time),4);
Endtime = Starttime + Interval;
globalstart = int32(round((Starttime - Flash_start_time)*10^(6),6));
globalend =  int32(globalstart + Interval*10^(6));
globalstart_buffer = globalstart - 100;
globalend_buffer = globalend + 100;
slidewindowmicroseconds = 100;% us
overlapmicroseconds = 50;  %us
SNR_Duke_matrix = zeros(int32((globalend - globalstart) / overlapmicroseconds - 1), 1,'single');
SNR_Hudson_matrix = zeros(int32((globalend - globalstart) / overlapmicroseconds - 1), 1,'single');
SNR_PS2_matrix = zeros(int32((globalend - globalstart) / overlapmicroseconds - 1), 1,'single');
SNR_PS3_matrix = zeros(int32((globalend - globalstart) / overlapmicroseconds - 1), 1,'single');
Abstime = zeros((globalend - globalstart) / overlapmicroseconds - 1, 1,'single');

%% Load data , Duke as base signal , the other three extra 200 data points extra to accomdate for shifting
%% Interpolate the data with 10 times time precision and High Pass filter it

sig_Duke_BGE_unfiltered = sm1.Duke_BGE';
sig_Duke_BGN_unfiltered = sm1.Duke_BGN';

sig_Hudson_BGE_unfiltered = sm1.Hudson_BGE';
sig_Hudson_BGN_unfiltered = sm1.Hudson_BGN';

sig_PS2_BGE_unfiltered = sm1.PS2_BGE';
sig_PS2_BGN_unfiltered = sm1.PS2_BGN';

sig_PS3_BGE_unfiltered = sm1.PS3_BGE';
sig_PS3_BGN_unfiltered = sm1.PS3_BGN';


%% Preflash standard deviation

Duke_BGE_preflashstd = single(std(sm1.Duke_BGE(1:1:int32(round((First_Pulse_time-Flash_start_time)*10^(6), 6)))));
Duke_BGN_preflashstd = single(std(sm1.Duke_BGN(1:1:int32(round((First_Pulse_time-Flash_start_time)*10^(6), 6)))));

Hudson_BGE_preflashstd = single(std(sm1.Hudson_BGE(1:1:int32(round((First_Pulse_time-Flash_start_time)*10^(6), 6)))));
Hudson_BGN_preflashstd = single(std(sm1.Hudson_BGN(1:1:int32(round((First_Pulse_time-Flash_start_time)*10^(6), 6)))));

PS2_BGE_preflashstd = single(std(sm1.PS2_BGE(1:1:int32(round((First_Pulse_time-Flash_start_time)*10^(6), 6)))));
PS2_BGN_preflashstd = single(std(sm1.PS2_BGN(1:1:int32(round((First_Pulse_time-Flash_start_time)*10^(6), 6)))));

PS3_BGE_preflashstd = single(std(sm1.PS3_BGE(1:1:int32(round((First_Pulse_time-Flash_start_time)*10^(6), 6)))));
PS3_BGN_preflashstd = single(std(sm1.PS3_BGN(1:1:int32(round((First_Pulse_time-Flash_start_time)*10^(6), 6)))));

%% Interpolation  

sig_Duke_BGE = single(interp(sig_Duke_BGE_unfiltered(globalstart:1:globalend),10));
sig_Duke_BGN = single(interp(sig_Duke_BGN_unfiltered(globalstart:1:globalend),10));

sig_Hudson_BGE = single(interp(sig_Hudson_BGE_unfiltered(globalstart:1:globalend),10));
sig_Hudson_BGN = single(interp(sig_Hudson_BGN_unfiltered(globalstart:1:globalend),10));

sig_PS2_BGE = single(interp(sig_PS2_BGE_unfiltered(globalstart:1:globalend),10));
sig_PS2_BGN = single(interp(sig_PS2_BGN_unfiltered(globalstart:1:globalend),10));

sig_PS3_BGE = single(interp(sig_PS3_BGE_unfiltered(globalstart:1:globalend),10));
sig_PS3_BGN = single(interp(sig_PS3_BGN_unfiltered(globalstart:1:globalend),10));


Totalpoints = int32((globalend - globalstart) / overlapmicroseconds) - 1;

Flashname = filename(1:1:19);
directoryname = sprintf('%s_Flash_%dus_slidingwindow_%dus_overlap_%9.6f-%9.6f_unfiltered_original_signal_only', Flashname,slidewindowmicroseconds, overlapmicroseconds, Starttime, Endtime);
mkdir(directoryname);
cd(directoryname);

parfor (i = 1:1:Totalpoints,20)
 [SNR_Duke,  SNR_Hudson,  SNR_PS2, SNR_PS3] = crosscorfunc((i - 1) * overlapmicroseconds*10 + 1, (i - 1) * overlapmicroseconds*10 + slidewindowmicroseconds*10 ,sig_Duke_BGE ,sig_Duke_BGN,sig_Hudson_BGE,sig_Hudson_BGN, sig_PS2_BGE, sig_PS2_BGN, sig_PS3_BGE,sig_PS3_BGN, First_Pulse_time, Duke_BGE_preflashstd, Duke_BGN_preflashstd, Hudson_BGE_preflashstd, Hudson_BGN_preflashstd, PS2_BGE_preflashstd, PS2_BGN_preflashstd, PS3_BGE_preflashstd, PS3_BGN_preflashstd);
 SNR_Duke_matrix(i) = SNR_Duke;
 SNR_Hudson_matrix(i) = SNR_Hudson;
 SNR_PS2_matrix(i) = SNR_PS2;
 SNR_PS3_matrix(i) = SNR_PS3;
 Abstime(i) = Starttime + single((single(i) - 1) * overlapmicroseconds*10^(-6));
 
end

filename = sprintf('%dus_slidingwindow_%dus_overlap_%9.6f-%9.6f.mat', slidewindowmicroseconds, overlapmicroseconds, Starttime, Endtime);
save(filename, 'SNR_Duke_matrix', 'SNR_Hudson_matrix', 'SNR_PS2_matrix', 'SNR_PS3_matrix', 'Abstime');
toc


function [SNR_Duke,  SNR_Hudson,  SNR_PS2, SNR_PS3] = crosscorfunc(startposition, Endposition, sig_Duke_BGE ,sig_Duke_BGN,sig_Hudson_BGE,sig_Hudson_BGN, sig_PS2_BGE, sig_PS2_BGN, sig_PS3_BGE,sig_PS3_BGN,Starttime, Duke_BGE_preflashstd, Duke_BGN_preflashstd, Hudson_BGE_preflashstd, Hudson_BGN_preflashstd, PS2_BGE_preflashstd, PS2_BGN_preflashstd, PS3_BGE_preflashstd, PS3_BGN_preflashstd)
    
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
        Duke_preflashstd = Duke_BGE_preflashstd;
        SNR_Duke = SNR_Duke_BGE;
    else
        sig_Duke  = sig_Duke_BGN(startposition:1:Endposition);
        Duke_preflashstd = Duke_BGN_preflashstd;
        SNR_Duke = SNR_Duke_BGN;
    end

    if (Hudson_signal == 1)
        sig_Hudson  = sig_Hudson_BGE(startposition:1:Endposition);
        Hudson_preflashstd = Hudson_BGE_preflashstd;
        SNR_Hudson = SNR_Hudson_BGE;
    else
        sig_Hudson  = sig_Hudson_BGN(startposition:1:Endposition);
        Hudson_preflashstd = Hudson_BGN_preflashstd;
        SNR_Hudson = SNR_Hudson_BGN;
    end

    if (PS2_signal == 1)
        sig_PS2  = sig_PS2_BGE(startposition:1:Endposition);
        PS2_preflashstd = PS2_BGE_preflashstd;
        SNR_PS2 = SNR_PS2_BGE;

    else
        sig_PS2  = sig_PS2_BGN(startposition:1:Endposition);
        PS2_preflashstd = PS2_BGN_preflashstd;
        SNR_PS2 = SNR_PS2_BGN;

    end

    if (PS3_signal == 1)
        sig_PS3  = sig_PS3_BGE(startposition:1:Endposition);
        PS3_preflashstd  = PS3_BGE_preflashstd;
        SNR_PS3 = SNR_PS3_BGE;

    else
        sig_PS3  = sig_PS3_BGN(startposition:1:Endposition);
        PS3_preflashstd  = PS3_BGN_preflashstd;
        SNR_PS3 = SNR_PS3_BGN;

    end

    [~, sm_max] = max(sig_Duke);
    [~, sm_min] = min(sig_Duke);
    [~, sn_max] = max(sig_Hudson);
    [~, sn_min] = min(sig_Hudson);
    [~, sl_max] = max(sig_PS2);
    [~, sl_min] = min(sig_PS2);
    [~, sp_max] = max(sig_PS3);
    [~, sp_min] = min(sig_PS3);

    if sm_max < sm_min

            % The polarity of first pulse of Duke is positive

            if sn_max < sn_min  %Hudson is positive
                sig_Hudson = 1*sig_Hudson;
            else
                sig_Hudson = -1*sig_Hudson;
            end

            if sl_max < sl_min
                sig_PS2 = 1*sig_PS2;
            else
                sig_PS2 = -1*sig_PS2;
            end

            if sp_max < sp_min
                sig_PS3 = 1*sig_PS3;
            else
                sig_PS3 = -1*sig_PS3;
            end

    else

            % The polarity of first pulse of Duke is negative

            if sn_max > sn_min  %Hudson is negative
                sig_Hudson = 1*sig_Hudson;
            else
                sig_Hudson = -1*sig_Hudson;
            end

            if sl_max > sl_min
                sig_PS2 = 1*sig_PS2;
            else
                sig_PS2 = -1*sig_PS2;
            end

            if sp_max > sp_min
                sig_PS3 = 1*sig_PS3;
            else
                sig_PS3 = -1*sig_PS3;
            end

    end

       % Max_std = [SNR_Duke,  SNR_Hudson,  SNR_PS2, SNR_PS3];
%         a = figure;
%         set(a, 'visible', 'off');
%         b = plot(1:1:length(sig_Duke),sig_Duke,'k-' , 1:1:length(sig_Duke), sig_Hudson, 'r', 1:1:length(sig_Duke), sig_PS2, 'b', 1:1:length(sig_Duke), sig_PS3, 'g');
%         legend('Duke', 'Hudson', 'Ps2', 'Ps3');
%         set(b(1), 'linewidth', 2);
%         set(b(2), 'linewidth', 0.75);
%         set(b(3), 'linewidth', 0.75);
%         set(b(4), 'linewidth', 0.75);
%         xlabel('time (0.1 us)');
%         ylabel('BGE/BGN Magnitude');
%         Abstimestart = Starttime + double(startposition-1)*0.0000001;
%         Abstimeend = Starttime + double(Endposition-1)*0.0000001;
%         
%         title(sprintf('Original signal characteristics %9.6fus to %9.6fus',Abstimestart, Abstimeend));
%         values = sprintf('%5.2f ', Max_std);
%         dim = [0.7 0.7 0.3 0.3];
%         annotation('textbox',dim , 'String', values);
%         name = sprintf('100us_slidingwindow_50us_overlap_%9.6fto_%9.6f.jpeg', Abstimestart, Abstimeend);
%         saveas(a,name);
%         delete(a);

    
end
