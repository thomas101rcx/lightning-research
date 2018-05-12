%% Cross-corr algorithm for lightning detection
tic
clear; format short e

mex fillgrid.c
mex timegrid.c
% mex linearto3d.c

XgridlimitUpper =   10;   %Km
XgridlimitLower =  -10;   %Km
YgridlimitUpper =   20;   %Km
YgridlimitLower =  -20;   %Km
Zgridlimit =        20;   %Km
xyGridsize   =     0.1;   %Km
zGridsize =        0.1;   %Km
xgridsize = length(XgridlimitLower:xyGridsize:XgridlimitUpper);
ygridsize = length(YgridlimitLower:xyGridsize:YgridlimitUpper);
zgridsize = length(0:zGridsize:Zgridlimit);
timegridsize = xgridsize*ygridsize*zgridsize;


%% Convert Long/Lat to relative x y coordinate with respect to Duke. Set Duke as base(0,0,0)

% Read in Site location text file

sitelocation = dlmread('./Fetched_LF/SiteLocations.txt', ' ', [0 1 3 2]);

NumStation = 4;

CenterStationLat =  35.97101;
CenterStationLon = -79.09433;

%Station1 Hudson
%Station2 Ps2
%Station3 Ps3

StationsLatLon = zeros(NumStation-1, 2, 'single');
Stations_xyz = zeros(NumStation-1, 3, 'single');

for i = 1:1:NumStation-1
    
    %% Manual enter coordinate
%     promptLat = sprintf('Enter Station %d Latitude\n', i);
%     promptLon = sprintf('Enter Station %d Longtitude\n', i);
%     answerLat = input(promptLat);
%     answerLon = input(promptLon);
    StationsLatLon(i,1) = sitelocation(i+1,1);
    StationsLatLon(i,2) = sitelocation(i+1,2);
    Dis_Center_Station_y = spheric_distance(CenterStationLat,CenterStationLon, StationsLatLon(i,1),CenterStationLon);
    Dis_Center_Station_x = spheric_distance(CenterStationLat,CenterStationLon, CenterStationLat, StationsLatLon(i,2));
    [StationLatPolairty, StationLonPolairty] = coorpolairty(CenterStationLat, CenterStationLon, StationsLatLon(i,1),StationsLatLon(i,2));
    
    Stations_xyz(i,1) = Dis_Center_Station_x.*StationLonPolairty;
    Stations_xyz(i,2) = Dis_Center_Station_y.*StationLatPolairty;
    Stations_xyz(i,3) = 0;
    
end


%% Calculate time difference with repsect to each grid

Dukex = 0;
Dukey = 0;
Dukez = 0;

syms x1 y1 x2 y2
eq1 = (1 - x1)*y1 == XgridlimitLower;
eq2 = (xgridsize - x1)*y1 == XgridlimitUpper;
sol1 = solve([eq1, eq2], [x1, y1]);
X1Solx = single(sol1.x1);
X1Soly = single(sol1.y1);

eq1 = (1 - x2)*y2 == YgridlimitLower;
eq2 = (ygridsize - x2)*y2 == YgridlimitUpper;
sol2 = solve([eq1, eq2], [x2, y2]);
Y1Solx = single(sol2.x2);
Y1Soly = single(sol2.y2);

BaseComb = nchoosek(NumStation,2);

Basex = [Dukex; Stations_xyz(:,1)]';

Basey = [Dukey; Stations_xyz(:,2)]';

Basez = [Dukez; Stations_xyz(:,3)]';

timedifference= zeros(xgridsize * ygridsize * zgridsize,1,'single');
timedifferencematrix = zeros(xgridsize * ygridsize * zgridsize, nchoosek(NumStation,2),'single');


parfor (j = 1:1:BaseComb,10)
    Combx = flipud(combnk(Basex,2));
    Comby = flipud(combnk(Basey,2));
    Combz = flipud(combnk(Basez,2));
    timedifferencematrix(:,j) = timegrid(timegridsize,xgridsize, ygridsize, X1Solx, X1Soly, Y1Solx, Y1Soly, Combx(j,1), Comby(j,1), Combz(j,1), Combx(j,2), Comby(j,2), Combz(j,2), timedifference);
end

ShortTimegrid1 = int16((-999:1:999)');
timedifferencematrix = int16(-1 * ((round(timedifferencematrix / (10 ^ (- 6)), 1) * 10)));

%TimeDiff 1 Duke - Hudson
%TimeDiff 2 Duke - Ps2e
%TimeDiff 3 Duke - Ps3
%TimeDiff 4 Hudson - Ps2
%TimeDiff 5 Hudson - Ps3
%TimeDiff 6 Ps2 - Ps3


%% Find starting point based on abs time
%% Load All necessary data


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
Interval = round(max(VHF_time)- min(VHF_time),4);
Endtime = Starttime + Interval;
globalstart = int32(round((Starttime - Flash_start_time)*10^(6),6));
globalend =  int32(globalstart + Interval*10^(6));

slidewindow = 100; %us
overlap = 25;  %us
Resultx = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');
Resulty = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');

Abstime = zeros((globalend - globalstart) / overlap - 1, 1,'single');
Image_maxoverstdxyz = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');
Image_max = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');
Image_std = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');
SNR_Duke_matrix = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');
SNR_Hudson_matrix = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');
SNR_PS2_matrix = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');
SNR_PS3_matrix = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');
Max_Duke_matrix = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');
Max_Hudson_matrix = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');
Max_PS2_matrix = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');
Max_PS3_matrix = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');
XCorr_Duke_Hudson_matrix = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');
XCorr_Duke_PS2_matrix = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');
XCorr_Duke_PS3_matrix = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');
XCorr_Hudson_PS2_matrix = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');
XCorr_Hudson_PS3_matrix = zeros(int32((globalend - globalstart) / overlap - 1),1,'single');
XCorr_PS2_PS3_matrix = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');
Max_Duke_index_matrix = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');
Min_Duke_index_matrix = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');
Max_Hudson_index_matrix = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');
Min_Hudson_index_matrix = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');
Max_PS2_index_matrix = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');
Min_PS2_index_matrix = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');
Max_PS3_index_matrix = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');
Min_PS3_index_matrix = zeros(int32((globalend - globalstart) / overlap - 1), 1,'single');


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

Totalpoints = int32((globalend - globalstart) / overlap) - 1;


%% Filtering for PS3
% remove noise at 150 KHz and 320 KHz
fs = 10^(7);
  
d = designfilt('bandstopiir','FilterOrder',6, ...
               'HalfPowerFrequency1',310000,'HalfPowerFrequency2',330000, ...
               'DesignMethod','butter','SampleRate',fs);
d1 = designfilt('bandstopiir','FilterOrder',6, ...
               'HalfPowerFrequency1',140000,'HalfPowerFrequency2',160000, ...
               'DesignMethod','butter','SampleRate',fs);
           
sig_PS3_BGE = filtfilt(d, double(sig_PS3_BGE));
sig_PS3_BGE = filtfilt(d1, double(sig_PS3_BGE));

sig_PS3_BGN = filtfilt(d, double(sig_PS3_BGN));
sig_PS3_BGN = filtfilt(d1, double(sig_PS3_BGN));


Flashname = filename(1:1:19);
directoryname = sprintf('%s_Flash_%dus_slidingwindow_%dus_overlap_%9.6f-%9.6f_All_window_Norm_Coeff_Abs', Flashname,slidewindow, overlap, Starttime, Endtime);
mkdir(directoryname);
cd(directoryname);
movefile ../fillgrid.mexmaci64

parfor (i = 1:1:Totalpoints-4,4)
    [X, Y, Z, Image_maxoverstd, Image_Maxvalue, STDImage, SNR_Duke, SNR_Hudson, SNR_PS2, SNR_PS3, Max_Duke, Max_Hudson, Max_PS2, Max_PS3, X_Corr_Duke_Hudson, X_Corr_Duke_PS2, X_Corr_Duke_PS3 ,X_Corr_Hudson_PS2, X_Corr_Hudson_PS3, X_Corr_PS2_PS3,Duke_max_index, Duke_min_index, Hudson_max_index, Hudson_min_index, PS3_max_index, PS3_min_index, PS2_max_index, PS2_min_index] = crosscorfunc((i - 1) * overlap*10 + 1, (i - 1) * overlap*10 + slidewindow*10 , ShortTimegrid1,sig_Duke_BGE ,sig_Duke_BGN,sig_Hudson_BGE,sig_Hudson_BGN, sig_PS2_BGE, sig_PS2_BGN, sig_PS3_BGE,sig_PS3_BGN ,xgridsize , ygridsize, zgridsize, timedifferencematrix, BaseComb, Duke_BGE_preflashstd, Duke_BGN_preflashstd, Hudson_BGE_preflashstd, Hudson_BGN_preflashstd, PS2_BGE_preflashstd, PS2_BGN_preflashstd, PS3_BGE_preflashstd, PS3_BGN_preflashstd, XgridlimitLower, XgridlimitUpper, YgridlimitLower, YgridlimitUpper, X1Solx, X1Soly, Y1Solx, Y1Soly, First_Pulse_time);
    Resultx(i) = (X(1) - X1Solx)*(X1Soly);
    Resulty(i) = (Y(1) - Y1Solx)*(Y1Soly);
    Resultz(i) = (Z(1) - 1) *0.1;
    Abstime(i) = Starttime + single((single(i) - 1) * overlap*10^(-6));
    Image_maxoverstdxyz(i) = Image_maxoverstd;
    Image_max(i) = Image_Maxvalue;
    Image_std(i) = STDImage;
    SNR_Duke_matrix(i) = SNR_Duke;
    SNR_Hudson_matrix(i) = SNR_Hudson;
    SNR_PS2_matrix(i) = SNR_PS2;
    SNR_PS3_matrix(i) = SNR_PS3;
    Max_Duke_matrix(i) = Max_Duke;
    Max_Hudson_matrix(i) = Max_Hudson;
    Max_PS2_matrix(i) = Max_PS2;
    Max_PS3_matrix(i) = Max_PS3;
    XCorr_Duke_Hudson_matrix(i) = X_Corr_Duke_Hudson;
    XCorr_Duke_PS2_matrix(i) = X_Corr_Duke_PS2;
    XCorr_Duke_PS3_matrix(i) = X_Corr_Duke_PS3;
    XCorr_Hudson_PS2_matrix(i) = X_Corr_Hudson_PS2;
    XCorr_Hudson_PS3_matrix(i) = X_Corr_Hudson_PS3;
    XCorr_PS2_PS3_matrix(i) = X_Corr_PS2_PS3;
    Max_Duke_index_matrix(i) = Duke_max_index;
    Min_Duke_index_matrix(i) = Duke_min_index;
    Max_Hudson_index_matrix(i) = Hudson_max_index;
    Min_Hudson_index_matrix(i) = Hudson_min_index;
    Max_PS2_index_matrix(i) = PS2_max_index;
    Min_PS2_index_matrix(i) = PS2_min_index;
    Max_PS3_index_matrix(i) = PS3_max_index;
    Min_PS3_index_matrix(i) = PS3_min_index;
end

filename = sprintf('%dus_slidingwindow_%dus_overlap_%9.6f-%9.6f_All_window.mat', slidewindow, overlap, Starttime, Endtime);
save(filename, 'Resultx', 'Resulty', 'Resultz', 'Image_maxoverstdxyz', 'Abstime', 'Image_max','Image_std', 'SNR_Duke_matrix', 'SNR_Hudson_matrix', 'SNR_PS2_matrix', 'SNR_PS3_matrix', 'Max_Duke_matrix', 'Max_Hudson_matrix', 'Max_PS2_matrix', 'Max_PS3_matrix', 'XCorr_Duke_Hudson_matrix', 'XCorr_Duke_PS2_matrix', 'XCorr_Duke_PS3_matrix', 'XCorr_Hudson_PS2_matrix', 'XCorr_Hudson_PS3_matrix', 'XCorr_PS2_PS3_matrix', 'Max_Duke_index_matrix', 'Min_Duke_index_matrix', 'Max_Hudson_index_matrix', 'Min_Hudson_index_matrix', 'Max_PS2_index_matrix', 'Min_PS2_index_matrix', 'Max_PS3_index_matrix', 'Min_PS3_index_matrix');
toc
function [X, Y, Z, Image_maxoverstd, Image_Maxvalue, STDImage, SNR_Duke, SNR_Hudson, SNR_PS2, SNR_PS3, Max_Duke, Max_Hudson, Max_PS2, Max_PS3, X_Corr_Duke_Hudson, X_Corr_Duke_PS2, X_Corr_Duke_PS3 ,X_Corr_Hudson_PS2, X_Corr_Hudson_PS3, X_Corr_PS2_PS3, ...
Duke_max_index, Duke_min_index, Hudson_max_index, Hudson_min_index, PS3_max_index, PS3_min_index, PS2_max_index, PS2_min_index] = crosscorfunc(startposition, Endposition, ShortTimegrid1, sig_Duke_BGE ,sig_Duke_BGN,sig_Hudson_BGE,sig_Hudson_BGN, sig_PS2_BGE, sig_PS2_BGN, sig_PS3_BGE,sig_PS3_BGN, xgridsize , ygridsize, zgridsize, timedifferencematrix, BaseComb,Duke_BGE_preflashstd, Duke_BGN_preflashstd, Hudson_BGE_preflashstd, Hudson_BGN_preflashstd, PS2_BGE_preflashstd, PS2_BGN_preflashstd, PS3_BGE_preflashstd, PS3_BGN_preflashstd,XgridlimitLower, XgridlimitUpper, YgridlimitLower, YgridlimitUpper, X1Solx, X1Soly, Y1Solx, Y1Soly, First_Pulse_time)




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
    
    
    % Polarity 
    
    % Find the max/min signal's index and return it 

%     if Duke_max_index < Duke_min_index
% 
%             % The polarity of first pulse of Duke is positive
% 
%             if Hudson_max_index < Hudson_min_index  %Hudson is positive
%                 sig_Hudson = 1*sig_Hudson;
%             else
%                 sig_Hudson = -1*sig_Hudson;
%             end
% 
%             if PS2_max_index < PS2_min_index
%                 sig_PS2 = 1*sig_PS2;
%             else
%                 sig_PS2 = -1*sig_PS2;
%             end
% 
%             if PS3_max_index < PS3_min_index
%                 sig_PS3 = 1*sig_PS3;
%             else
%                 sig_PS3 = -1*sig_PS3;
%             end
% 
%     else
% 
%             % The polarity of first pulse of Duke is negative
% 
%             if Hudson_max_index > Hudson_min_index  %Hudson is negative
%                 sig_Hudson = 1*sig_Hudson;
%             else
%                 sig_Hudson = -1*sig_Hudson;
%             end
% 
%             if PS2_max_index > PS2_min_index
%                 sig_PS2 = 1*sig_PS2;
%             else
%                 sig_PS2 = -1*sig_PS2;
%             end
% 
%             if PS3_max_index > PS3_min_index
%                 sig_PS3 = 1*sig_PS3;
%             else
%                 sig_PS3 = -1*sig_PS3;
%             end
% 
%     end
    
   

  
%     if SNR_Duke > 10 && SNR_Hudson > 10 && SNR_PS2 > 10 && SNR_PS3 > 10
         %%% Long term goal : Make .mex file for the following
         
       
%  
%         f1 = figure(3);
%         subplot(2,2,[1 2]);
%         set(f1, 'visible', 'off');
%         b = plot(1:1:length(sig_Duke),sig_Duke,'k-' , 1:1:length(sig_Duke), sig_Hudson, 'r', 1:1:length(sig_Duke), sig_PS2, 'b', 1:1:length(sig_Duke), sig_PS3, 'g');
%         legend('Duke', 'Hudson', 'Ps2', 'Ps3');
%         set(b(1), 'linewidth', 2);
%         set(b(2), 'linewidth', 0.75);
%         set(b(3), 'linewidth', 0.75);
%         set(b(4), 'linewidth', 0.75);
%         xlabel('time (0.1 us)');
%         ylabel('BGE/BGN Magnitude');
%         Abstimestart = First_Pulse_time + double(startposition-1)*0.0000001;
%         Abstimeend = First_Pulse_time + double(Endposition-1)*0.0000001;
%         title(sprintf('Original signal characteristics %9.6fus to %9.6fus',Abstimestart, Abstimeend), 'Fontsize', 9);
%         Max_std = [SNR_Duke,  SNR_Hudson,  SNR_PS2, SNR_PS3];
%         values = sprintf('%5.2f ', Max_std);
%         dim = [0.7 0.7 0.3 0.3];
%         annotation('textbox',dim , 'String', values);

        %% Generate Ixy matrix
        
        
   
        Ixymatrix = zeros(length(ShortTimegrid1), BaseComb,'single');
        
        
        [Ixymatrix(:,1), ~] = xcorr(sig_Hudson, sig_Duke, 'coeff');
        [Ixymatrix(:,2), ~] = xcorr(sig_PS2, sig_Duke, 'coeff');
        [Ixymatrix(:,3), ~] = xcorr(sig_PS3, sig_Duke, 'coeff');
        [Ixymatrix(:,4), ~] = xcorr(sig_PS2, sig_Hudson, 'coeff');
        [Ixymatrix(:,5), ~] = xcorr(sig_PS3, sig_Hudson, 'coeff');
        [Ixymatrix(:,6), ~] = xcorr(sig_PS3, sig_PS2, 'coeff');
        
        X_Corr_Duke_Hudson = max(abs(Ixymatrix(:,1)));
        X_Corr_Duke_PS2 = max(abs(Ixymatrix(:,2)));
        X_Corr_Duke_PS3 = max(abs(Ixymatrix(:,3)));
        X_Corr_Hudson_PS2 = max(abs(Ixymatrix(:,4)));
        X_Corr_Hudson_PS3 = max(abs(Ixymatrix(:,5)));
        X_Corr_PS2_PS3 = max(abs(Ixymatrix(:,6)));

        
%         subplot(2,2,3);
%         c = plot(1:1:length(Ixymatrix(:,1)),Ixymatrix(:,1),'k-' , 1:1:length(Ixymatrix(:,1)), Ixymatrix(:,2)+2, 'r', 1:1:length(Ixymatrix(:,1)), Ixymatrix(:,3) -2 , 'b', 1:1:length(Ixymatrix(:,1)), Ixymatrix(:,4) +4, 'g',1:1:length(Ixymatrix(:,1)), Ixymatrix(:,5) -4, 'y', 1:1:length(Ixymatrix(:,1)), Ixymatrix(:,6) + 6 ,'m');
%         set(c(1), 'linewidth', 0.75);
%         set(c(2), 'linewidth', 0.75);
%         set(c(3), 'linewidth', 0.75);
%         set(c(4), 'linewidth', 0.75);
%         set(c(5), 'linewidth', 0.75);
%         set(c(6), 'linewidth', 0.75);
%         xlabel('seconds', 'FontSize', 13);
%         ylabel('Magnitude of cross-correlating', 'FontSize', 13);
%         lgd = legend('D/H', 'D/PS2', 'D/PS3', 'H/PS2', 'H/PS3', 'PS3/PS2');
%         set(lgd, 'Fontsize', 7);
%         set(lgd, 'Position', [0.8 0 0.2 0.2]); 
%         title('Cross-correlation distribution over 100 \museconds', 'FontSize', 9);


        %% Fill in the grids
        SumIxyz = zeros(xgridsize*ygridsize*zgridsize,1,'single');

        parfor (i = 1:1:BaseComb, 4)
            Ixyzmex = fillgrid(timedifferencematrix(:,i), abs(Ixymatrix(:,i)));
            SumIxyz =  SumIxyz+ Ixyzmex ;
        end


        [Image_Maxvalue, Index] = max(SumIxyz(:));
        [X, Y, Z] = linearto3d(Index, xgridsize, ygridsize);
        STDImage = std(SumIxyz(:));
        Image_maxoverstd = Image_Maxvalue/STDImage;
       
        
         
        E = reshape(SumIxyz, [xgridsize ygridsize zgridsize]);
        
        x = XgridlimitLower:0.1:XgridlimitUpper;
        y = YgridlimitLower:0.1:YgridlimitUpper;
        
        X1 = (X(1) - X1Solx)*(X1Soly);
        Y1 = (Y(1) - Y1Solx)*(Y1Soly);
        Z1 = (Z(1) - 1) *0.1;
%         
%         subplot(2,2,4);
%         hold on
%         imagesc(x,y,E(:,:,Z)');
%         plot(X1,Y1, 'k*', 'markersize', 10);
%         set(gca,'YDir','normal');
%         xlabel('X Direction (km)','FontSize', 13);
%         ylabel('Y Direction (km)','FontSize', 13);
%         Frametime = double(startposition)*0.0000001 + double(First_Pulse_time);
%         title(sprintf('X:%g km, Y:%g km, Z:%g Km,%f time', X1, Y1, Z1, Frametime), 'FontSize', 9);
%         xlim([XgridlimitLower XgridlimitUpper]);
%         ylim([YgridlimitLower YgridlimitUpper]);
%         xticks(linspace(XgridlimitLower, XgridlimitUpper, 5))
%         yticks(linspace(YgridlimitLower, YgridlimitUpper, 5))
%         colorbar
%         hold off
%         saveas(f1,sprintf('%f time.jpeg', Frametime), 'jpeg');
%         delete(f1);
%         clf
        
%     else
%         X = -10000;
%         Y = -10000;
%         Z = -10000;
%         maxoverstd = -10000;
%     end
    
end

%% Calculate distance between two coordinates
function Dists = spheric_distance(lat1, lon1, lat2, lon2)
    %%%
    %%% This function calculate the spheric distance between two location on
    %%% earth.   kilometers
    %%%
    lat1 = double(lat1);
    lon1 = double(lon1);
    lat2 = double(lat2);
    lon2 = double(lon2);

    R = 6371.0;
    pdeg2rad = pi / 180.0;

    lat1 = lat1 * pdeg2rad;
    lat2 = lat2 * pdeg2rad;
    lon1 = lon1 * pdeg2rad;
    lon2 = lon2 * pdeg2rad;

    Dists = (R * acos(sin(lat1) .* sin(lat2) + cos(lat1) .* cos(lat2) .* cos(lon2 - lon1)));

end

%% Convert 1D index to 3D index

function [i ,j, k] = linearto3d(index, xgridsize, ygridsize)

    i = mod(index-1,xgridsize)+1;
    j = mod((index - i)/xgridsize, ygridsize)+1;
    k = ((index - i)/xgridsize-j+1)/ygridsize + 1;

end

%% Determine other stations' relative location to the center station

function [StationLatPolairty, StationLonPolairty] = coorpolairty(CenterStationLat, CenterStationLon, StationLat, StationLon)
    
    if StationLat- CenterStationLat > 0
        StationLatPolairty = 1;
    elseif StationLat - CenterStationLat < 0
        StationLatPolairty = - 1;
    else
        StationLatPolairty = 1;
    end
    
    if StationLon - CenterStationLon > 0
        StationLonPolairty = 1;
    elseif StationLon- CenterStationLon < 0 
        StationLonPolairty = -1;
    else
        StationLonPolairty = 1;
    end

end

