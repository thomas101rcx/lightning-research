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

ShortTimegrid1 = int16((-1000:1:1000)');
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

Starttime = First_Pulse_time;
Interval = 0.05
%round(max(VHF_time)- min(VHF_time),4);
Endtime = Starttime + Interval;
globalstart = int32(round((Starttime - Flash_start_time)*10^(6),6));
globalend =  int32(globalstart + Interval*10^(6));
globalstart_buffer = globalstart - 100;
globalend_buffer = globalend + 100;
slidewindow = 50; %us
slidewindowmicroseconds = slidewindow;
overlap = 25;  %us
overlapmicroseconds = overlap;
Resultx = zeros(int32((globalend - globalstart) / overlapmicroseconds - 1), 1,'single');
Resulty = zeros(int32((globalend - globalstart) / overlapmicroseconds - 1), 1,'single');
Resultz = zeros(int32((globalend - globalstart) / overlapmicroseconds - 1), 1,'single');
Abstime = zeros((globalend - globalstart) / overlapmicroseconds - 1, 1,'single');
maxoverstdxyz = zeros(int32((globalend - globalstart) / overlapmicroseconds - 1), 1,'single');

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

sig_Duke_BGE = single(interp(sig_Duke_BGE_unfiltered(globalstart_buffer:1:globalend_buffer),10));
sig_Duke_BGN = single(interp(sig_Duke_BGN_unfiltered(globalstart_buffer:1:globalend_buffer),10));

sig_Hudson_BGE = single(interp(sig_Hudson_BGE_unfiltered(globalstart_buffer:1:globalend_buffer),10));
sig_Hudson_BGN = single(interp(sig_Hudson_BGN_unfiltered(globalstart_buffer:1:globalend_buffer),10));

sig_PS2_BGE = single(interp(sig_PS2_BGE_unfiltered(globalstart_buffer:1:globalend_buffer),10));
sig_PS2_BGN = single(interp(sig_PS2_BGN_unfiltered(globalstart_buffer:1:globalend_buffer),10));

sig_PS3_BGE = single(interp(sig_PS3_BGE_unfiltered(globalstart_buffer:1:globalend_buffer),10));
sig_PS3_BGN = single(interp(sig_PS3_BGN_unfiltered(globalstart_buffer:1:globalend_buffer),10));


Totalpoints = int32((globalend - globalstart) / overlapmicroseconds) - 1;

Flashname = filename(1:1:19);
directoryname = sprintf('%s_Flash_%dus_slidingwindow_%dus_overlap_%9.6f-%9.6f_Interp', Flashname,slidewindowmicroseconds, overlapmicroseconds, Starttime, Endtime);
mkdir(directoryname);
cd(directoryname);
movefile ../fillgrid.mexmaci64


parfor (i = 1:1:Totalpoints, 4)
    [X, Y, Z, maxoverstd] = crosscorfunc((i - 1) * overlap*10 + 1, (i - 1) * overlap*10 + slidewindow*10 , ShortTimegrid1,sig_Duke_BGE ,sig_Duke_BGN,sig_Hudson_BGE,sig_Hudson_BGN, sig_PS2_BGE, sig_PS2_BGN, sig_PS3_BGE,sig_PS3_BGN ,xgridsize , ygridsize, zgridsize, timedifferencematrix, BaseComb, Duke_BGE_preflashstd, Duke_BGN_preflashstd, Hudson_BGE_preflashstd, Hudson_BGN_preflashstd, PS2_BGE_preflashstd, PS2_BGN_preflashstd, PS3_BGE_preflashstd, PS3_BGN_preflashstd, XgridlimitLower, XgridlimitUpper, YgridlimitLower, YgridlimitUpper, X1Solx, X1Soly, Y1Solx, Y1Soly, First_Pulse_time);
    Resultx(i) = (X(1) - X1Solx)*(X1Soly);
    Resulty(i) = (Y(1) - Y1Solx)*(Y1Soly);
    Resultz(i) = (Z(1) - 1) *0.1;
    Abstime(i) = Starttime + single((single(i) - 1) * overlapmicroseconds*10^(-6));
    maxoverstdxyz(i) = maxoverstd;
end

filename = sprintf('%dus_slidingwindow_%dus_overlap_%9.6f-%9.6f.mat', slidewindowmicroseconds, overlapmicroseconds, Starttime, Endtime);
save(filename, 'Resultx', 'Resulty', 'Resultz', 'maxoverstdxyz', 'Abstime');
toc
function [X, Y, Z,maxoverstd] = crosscorfunc(startposition, Endposition, ShortTimegrid1, sig_Duke_BGE ,sig_Duke_BGN,sig_Hudson_BGE,sig_Hudson_BGN, sig_PS2_BGE, sig_PS2_BGN, sig_PS3_BGE,sig_PS3_BGN, xgridsize , ygridsize, zgridsize, timedifferencematrix, BaseComb,Duke_BGE_preflashstd, Duke_BGN_preflashstd, Hudson_BGE_preflashstd, Hudson_BGN_preflashstd, PS2_BGE_preflashstd, PS2_BGN_preflashstd, PS3_BGE_preflashstd, PS3_BGN_preflashstd,XgridlimitLower, XgridlimitUpper, YgridlimitLower, YgridlimitUpper, X1Solx, X1Soly, Y1Solx, Y1Soly, First_Pulse_time)

    buffersize = 2000;
     
    
%     Determine which signal to use BGE or BGN based on the max(abs) value of the signal
%
%     Choose Duke as signal polarity reference
%
%     Signal 1 = BGE;
%     Signal 2 = BGN;



    SNR_Duke_BGE = sig_Duke_BGE/Duke_BGE_preflashstd;
    SNR_Duke_BGN = sig_Duke_BGN/Duke_BGN_preflashstd;
    
    SNR_Hudson_BGE = sig_Hudson_BGE/Hudson_BGE_preflashstd;
    SNR_Hudson_BGN = sig_Hudson_BGN/Hudson_BGN_preflashstd;
    
    SNR_PS2_BGE = sig_PS2_BGE/PS2_BGE_preflashstd;
    SNR_PS2_BGN = sig_PS2_BGN/PS2_BGN_preflashstd;

    SNR_PS3_BGE = sig_PS3_BGE/PS3_BGE_preflashstd;
    SNR_PS3_BGN = sig_PS3_BGN/PS3_BGN_preflashstd;
    
   
    Duke_signal_matrix_max = [max(abs(SNR_Duke_BGE(startposition:1:Endposition))), max(abs(SNR_Duke_BGN(startposition:1:Endposition)))];
    Hudson_signal_matrix_max = [max(abs(SNR_Hudson_BGE(startposition:1:Endposition+buffersize))),max(abs(SNR_Hudson_BGN(startposition:1:Endposition+buffersize)))];
    PS2_signal_matrix_max = [max(abs(SNR_PS2_BGE(startposition:1:Endposition+buffersize))),max(abs(SNR_PS2_BGN(startposition:1:Endposition+buffersize)))];
    PS3_signal_matrix_max = [max(abs(SNR_PS3_BGE(startposition:1:Endposition+buffersize))),max(abs(SNR_PS3_BGN(startposition:1:Endposition+buffersize)))];

    [~, Duke_signal] = max(Duke_signal_matrix_max);
    [~, Hudson_signal] = max(Hudson_signal_matrix_max);
    [~, PS2_signal] = max(PS2_signal_matrix_max);
    [~, PS3_signal] = max(PS3_signal_matrix_max);
    
   if (Duke_signal == 1)
        sig_Duke  = sig_Duke_BGE(startposition:1:Endposition);
        SNR_Duke = SNR_Duke_BGE;
    else
        sig_Duke  = sig_Duke_BGN(startposition:1:Endposition);
        SNR_Duke = SNR_Duke_BGN;
    end

    if (Hudson_signal == 1)
        sig_Hudson  = sig_Hudson_BGE(startposition:1:Endposition);
        SNR_Hudson = SNR_Hudson_BGE;
    else
        sig_Hudson  = sig_Hudson_BGN(startposition:1:Endposition);
        SNR_Hudson = SNR_Hudson_BGN;
    end

    if (PS2_signal == 1)
        sig_PS2  = sig_PS2_BGE(startposition:1:Endposition);
        SNR_PS2 = SNR_PS2_BGE;

    else
        sig_PS2  = sig_PS2_BGN(startposition:1:Endposition);
        SNR_PS2 = SNR_PS2_BGN;

    end

    if (PS3_signal == 1)
        sig_PS3  = sig_PS3_BGE(startposition:1:Endposition);
        SNR_PS3 = SNR_PS3_BGE;

    else
        sig_PS3  = sig_PS3_BGN(startposition:1:Endposition);
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

    %Window the signal so that we get the best resolution
  %  vq_sm = single(vq_sm.* gausswin(length(vq_sm)));
    
    % Determine if this window is valid or not by looking at the average of
    % the abs value of that window , if it's greater than the backgroud
    % noise by 1.2 times, count as valid. if not, then it's not valid
    
    % Mean/Average method 
   % vq_sm_sig = mean(abs(vq_sm));
    
    % Standard deviation method
    
%     vq_sn_sig_std = std(vq_sn);
%     vq_sl_sig_std = std(vq_sl);
%     vq_sp_sig_std = std(vq_sp);
    
    if  SNR_Duke > 10 && SNR_Hudson > 10 && SNR_PS2 > 10 && SNR_PS3 > 10
         %%% Long term goal : Make .mex file for the following

        %% Generate Ixy matrix
   
        Ixymatrix = zeros(length(ShortTimegrid1), BaseComb,'single');

        %% Dot product calculation based on shiftlength:
        for y = 1:1:length(ShortTimegrid1)
            startpoint = 1001;
            endpoint = 1500;

            sn_test_shift = sig_Hudson(startpoint + ShortTimegrid1(y):1:endpoint + ShortTimegrid1(y),1); %Hudson
            sl_test_shift = sig_PS2(startpoint + ShortTimegrid1(y):1:endpoint + ShortTimegrid1(y),1); %Ps2
            sp_test_shift = sig_PS3(startpoint + ShortTimegrid1(y):1:endpoint + ShortTimegrid1(y),1); %Ps3

            Ixymatrix(y,1) = sig_Duke'*sn_test_shift;   %Duke Hudson
            Ixymatrix(y,2) = sig_Duke'*sl_test_shift;   %Duke Ps2
            Ixymatrix(y,3) = sig_Duke'*sp_test_shift;   %Duke Ps3
            Ixymatrix(y,4) = sig_Hudson(startpoint:1:endpoint)'*sl_test_shift; %Hudson Ps2
            Ixymatrix(y,5) = sig_Hudson(startpoint:1:endpoint)'*sp_test_shift; % Hudson Ps3
            Ixymatrix(y,6) = sig_PS2(startpoint:1:endpoint)'*sp_test_shift; % Ps2 Ps3

        end
        
        plot(Ixymatrix(:,1));
        plot(Ixymatrix(:,2));
        plot(Ixymatrix(:,3));
        plot(Ixymatrix(:,4));

        %% Fill in the grids
        SumIxyz = zeros(xgridsize*ygridsize*zgridsize,1,'single');

        % Remember to change fillgrid.c to 2001
        parfor (i = 1:1:BaseComb,4)
            Ixyzmex = fillgrid(timedifferencematrix(:,i), Ixymatrix(:,i));
            SumIxyz =  SumIxyz+ Ixyzmex/max(Ixyzmex) ;
        end


        [Maxvalue, Index] = max(SumIxyz(:));
        [X, Y, Z] = linearto3d(Index, xgridsize, ygridsize);
        STDImage = std(SumIxyz(:));
        maxoverstd = Maxvalue/STDImage;
        E = reshape(SumIxyz, [xgridsize ygridsize zgridsize]);
        
        x = XgridlimitLower:0.1:XgridlimitUpper;
        y = YgridlimitLower:0.1:YgridlimitUpper;
        
        X1 = (X(1) - X1Solx)*(X1Soly);
        Y1 = (Y(1) - Y1Solx)*(Y1Soly);
        Z1 = (Z(1) - 1) *0.1;
        
        f10 = figure(10);
        hold on
        set(f10, 'Visible', 'off')
        imagesc(x,y,E(:,:,Z)');
        plot(X1,Y1, 'k*', 'markersize', 10);
        set(gca,'YDir','normal');
        xlabel('X Direction (km)','FontSize', 17);
        ylabel('Y Direction (km)','FontSize', 17);
        Frametime = double(startposition)*0.0000001 + double(First_Pulse_time);
        title(sprintf('Horizontal Cross-section at %g km %f time', Z1, Frametime), 'FontSize', 20);
        xlim([XgridlimitLower XgridlimitUpper]);
        ylim([YgridlimitLower YgridlimitUpper]);
        xticks(linspace(XgridlimitLower, XgridlimitUpper, 5))
        yticks(linspace(YgridlimitLower, YgridlimitUpper, 5))
        colorbar
        hold off
        saveas(f10,sprintf('All-Stations-Horizontal Cross-section %f time at %g km.jpeg', Z1, Frametime), 'jpeg');
        delete(f10);
        clf
        
    else
        X = -10000;
        Y = -10000;
        Z = -10000;
        maxoverstd = -10000;
    end
    
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

