%% Cross-corr algorithm for lightning detection
clear; format short e
tic;

XgridlimitUpper =  5;    %Km
XgridlimitLower = -10;   %Km
YgridlimitUpper =  15;   %Km
YgridlimitLower =   0;   %Km
Zgridlimit =       20;   %Km
xyGridsize   =    0.1;   %Km
zGridsize =       0.2;   %Km
xgridsize = length(XgridlimitLower:xyGridsize:XgridlimitUpper);
ygridsize = length(YgridlimitLower:xyGridsize:YgridlimitUpper);
zgridsize = length(0:zGridsize:Zgridlimit);

timedifference1 = zeros(xgridsize * ygridsize * zgridsize, 1);
timedifference2 = zeros(xgridsize * ygridsize * zgridsize, 1);
timedifference3 = zeros(xgridsize * ygridsize * zgridsize, 1);
timedifference4 = zeros(xgridsize * ygridsize * zgridsize, 1);
timedifference5 = zeros(xgridsize * ygridsize * zgridsize, 1);
timedifference6 = zeros(xgridsize * ygridsize * zgridsize, 1);


%% Convert Long/Lat to relative x y coordinate with respect to Duke. Set Duke as base(0,0,0)

DukeLat = 35.97101;
DukeLon = -79.09433;

HudsonLat =  36.00423;
HudsonLon = -78.94033;

Ps2Lat = 35.943668;
Ps2Lon = -79.321576;

Ps3Lat = 35.807490;
Ps3Lon = -78.917833;


% Duke to Hudson
Dis_Duke_Hudson   = spheric_distance(DukeLat,DukeLon, HudsonLat, HudsonLon);
Dis_Duke_Hudson_y = spheric_distance(DukeLat,DukeLon, HudsonLat,DukeLon);
Dis_Duke_Hudson_x = spheric_distance(DukeLat,DukeLon, DukeLat, HudsonLon);

Hudson_xyz = [Dis_Duke_Hudson_x, Dis_Duke_Hudson_y, 0];

% Duke to Ps2

Dis_Duke_Ps2   = spheric_distance(DukeLat,DukeLon, Ps2Lat,   Ps2Lon);
Dis_Duke_Ps2_y = spheric_distance(DukeLat,DukeLon, Ps2Lat, DukeLon);
Dis_Duke_Ps2_x = spheric_distance(DukeLat,DukeLon, DukeLat, Ps2Lon);

Ps2_xyz = [- Dis_Duke_Ps2_x, -Dis_Duke_Ps2_y, 0];

% Duke to PS3

Dis_Duke_Ps3   = spheric_distance(DukeLat,DukeLon, Ps3Lat,  Ps3Lon);
Dis_Duke_Ps3_y = spheric_distance(DukeLat,DukeLon, Ps3Lat, DukeLon);
Dis_Duke_Ps3_x = spheric_distance(DukeLat,DukeLon, DukeLat, Ps3Lon);

Ps3_xyz = [Dis_Duke_Ps3_x, -Dis_Duke_Ps3_y, 0];

Hudsonx = Hudson_xyz(1);
Hudsony = Hudson_xyz(2);
Hudsonz = Hudson_xyz(3);

Ps2x = Ps2_xyz(1);
Ps2y = Ps2_xyz(2);
Ps2z = Ps2_xyz(3);

Ps3x = Ps3_xyz(1);
Ps3y = Ps3_xyz(2);
Ps3z = Ps3_xyz(3);


syms x1 y1 x2 y2  
eq1 = (1 - x1)*y1 == XgridlimitLower;
eq2 = (xgridsize - x1)*y1 == XgridlimitUpper;
sol1 = solve([eq1, eq2], [x1, y1]);
X1Solx = double(sol1.x1);
X1Soly = double(sol1.y1);

eq1 = (1 - x2)*y2 == YgridlimitLower;
eq2 = (ygridsize - x2)*y2 == YgridlimitUpper;
sol2 = solve([eq1, eq2], [x2, y2]);
Y1Solx = double(sol2.x2);
Y1Soly = double(sol2.y2);


%% Calculate time difference with repsect to each grid

parfor i = 1:1: xgridsize* ygridsize* zgridsize
    [X, Y, Z] = ind2sub([xgridsize ygridsize zgridsize], i);
    X1 = (X - ceil(xygridsize/2))*(XYgridlimit/floor(xygridsize/2));
    Y1 = (Y - ceil(xygridsize/2))*(XYgridlimit/floor(xygridsize/2));
    Z1 = (Z - 1)*(Gridsize);
    % Grid to Duke - Grid to other stations
    timedifference1(i) = calcudis(0, 0, 0, Hudsonx, Hudsony, 0, X1, Y1, Z1);
    timedifference2(i) = calcudis(0, 0, 0, Ps2x, Ps2y, 0, X1, Y1, Z1);
    timedifference3(i) = calcudis(0, 0, 0, Ps3x, Ps3y, 0, X1, Y1, Z1);
  
    % Grid to Hudson - Grid to other stations 
    
    timedifference4(i) = calcudis(Hudsonx, Hudsony, 0, Ps2x, Ps2y, 0, X1, Y1, Z1);
    timedifference5(i) = calcudis(Hudsonx, Hudsony, 0, Ps3x, Ps3y, 0, X1, Y1, Z1);
  
    % Grid to PS2 - Grid to other stations 
    
    timedifference6(i) = calcudis(Ps2x, Ps2y, 0, Ps3x, Ps3y, 0, X1, Y1, Z1);
end

ShortTimegrid1 = (-200:1:200);

TimedifferenceGrid1 =  -1 * ((round(timedifference1(:) / (10 ^ (- 6)), 0) * 1));
TimedifferenceGrid2 =  -1 * ((round(timedifference2(:) / (10 ^ (- 6)), 0) * 1));
TimedifferenceGrid3 =  -1 * ((round(timedifference3(:) / (10 ^ (- 6)), 0) * 1));
TimedifferenceGrid4 =  -1 * ((round(timedifference4(:) / (10 ^ (- 6)), 0) * 1));
TimedifferenceGrid5 =  -1 * ((round(timedifference5(:) / (10 ^ (- 6)), 0) * 1));
TimedifferenceGrid6 =  -1 * ((round(timedifference6(:) / (10 ^ (- 6)), 0) * 1));

%TimeDiff 1 Duke -Hudson 
%TimeDiff 2 Duke -Ps2 
%TimeDiff 3 Duke -Ps3
%TimeDiff 4 Hudson -Ps2
%TimeDiff 5 Hudson -Ps3 
%TimeDiff 6 Ps2 - Ps3 




globalstart = 196300;
globalend =   296300;
globalstart_buffer = globalstart - 200;
globalend_buffer = globalend + 200;
slidewindow = 50; %us
slidewindowmicroseconds = slidewindow;
overlap = 25; %  us
overlapmicroseconds = overlap;
Totaltime = (globalend-globalstart )/ 1000; 
Resultx = zeros((globalend - globalstart) / overlapmicroseconds - 1, 1);
Resulty = zeros((globalend - globalstart) / overlapmicroseconds - 1, 1);
Resultz = zeros((globalend - globalstart) / overlapmicroseconds - 1, 1);
Abstime = zeros((globalend - globalstart) / overlapmicroseconds - 1, 1);
maxoverstdxyz = zeros((globalend - globalstart) / overlapmicroseconds - 1, 1);


%% Load data , Duke as base signal , the other three extra 200 data points extra to accomdate for shifting
sm1 = load('./Fetched_LF/2017-07-08_20-58-10_Flash.mat');


fc = 100000;
fs = 10^(6);
order = 10;
%High Pass filtering with butterworth response of order 10 
[b1,a1] = butter(order, fc/(fs/2),'high');

sm_BGE_filtered = filtfilt(b1,a1,sm1.Duke_BGE);
sn_BGE_filtered = filtfilt(b1,a1,sm1.Hudson_BGE);
sl_BGE_filtered = filtfilt(b1,a1,sm1.PS2_BGE);
sp_BGE_filtered = filtfilt(b1,a1,sm1.PS3_BGE);


sm_BGN_filtered = filtfilt(b1,a1,sm1.Duke_BGN);
sn_BGN_filtered = filtfilt(b1,a1,sm1.Hudson_BGN);
sl_BGN_filtered = filtfilt(b1,a1,sm1.PS2_BGN);
sp_BGN_filtered = filtfilt(b1,a1,sm1.PS3_BGN);


sm_BGE = sm_BGE_filtered(globalstart:1:globalend);
sm_BGN = sm_BGN_filtered(globalstart:1:globalend);

sn_BGE = sn_BGE_filtered(globalstart_buffer:1:globalend_buffer);
sn_BGN = sn_BGN_filtered(globalstart_buffer:1:globalend_buffer);

sl_BGE = sl_BGE_filtered(globalstart_buffer:1:globalend_buffer);
sl_BGN = sl_BGN_filtered(globalstart_buffer:1:globalend_buffer);

sp_BGE = sp_BGE_filtered(globalstart_buffer:1:globalend_buffer);
sp_BGN = sp_BGN_filtered(globalstart_buffer:1:globalend_buffer);


% t = 0:1/fs:2-1/fs;
% X = fft(sm1.PS2_BGE);
% 
% n = length(sm1.PS2_BGE); % length(x) gives the array length of signal x
% 
% c = (-1 * fs) / 2:fs / n:fs / 2 - fs / n; % It generates the frequency series to plot X in frequency domain
% 
% figure(2);
% 
% subplot(4, 1, 1),plot(t,sm1.PS2_BGE); % This subplot shows the signal x vs. time series t
% subplot(4, 1 ,2),plot(c,fftshift(abs(X))); % This subplot shows the Fourier spectrum of x with zero frequency component shifted to center
% subplot(4, 1, 3),plot(c,angle(X)); % This subplot shows the phase distribution of X (Fourier transform of x)
% subplot(4,1 ,4),plot(c,real(X)); % This subplot shows the real component of X spectrum
% plot(f,mx);
% title('Power Spectrum of a Sine Wave'); xlabel('Frequency (Hz)'); ylabel('Power');
% figure(1)
% plot(sm_BGE)
% hold on
% plot(sm_BGE_filtered)
% hold off
% figure(2)
% plot(sm_BGN)
% hold on
% plot(sm_BGN_filtered)
% hold off
% 
% figure(3)
% plot(sn_BGE)
% hold on
% plot(sn_BGE_filtered)
% hold off
% figure(4)
% plot(sn_BGN)
% hold on
% plot(sn_BGN_filtered)
% hold off
% 
% 
% 
% figure(5)
% plot(sl_BGE)
% hold on
% plot(sl_BGE_filtered)
% hold off
% figure(6)
% plot(sl_BGN)
% hold on
% plot(sl_BGN_filtered)
% hold off
% 
% 
% 
% figure(7)
% plot(sp_BGE)
% hold on
% plot(sp_BGE_filtered)
% hold off
% figure(8)
% plot(sp_BGN)
% hold on
% plot(sp_BGN_filtered)
% hold off
% 
% 
% 
% figure(9)
% plot(sq_BGE)
% hold on
% plot(sq_BGE_filtered)
% hold off
% figure(10)
% plot(sq_BGN)
% hold on
% plot(sq_BGN_filtered)
% hold off

directoryname = sprintf('2017-07-08_20-58-10_Flash_%dus_slidingwindow_%dus_overlap_%dms_data_Image_method_HighPass', slidewindowmicroseconds, overlapmicroseconds, Totaltime);
mkdir(directoryname);
cd(directoryname);
Totalpoints = ((globalend - globalstart) / overlapmicroseconds) - 1;

parfor i = 1:1:Totalpoints
    [X, Y, Z, maxoverstd] = crosscorfunc((i - 1) * overlap + 1, (i - 1) * overlap + slidewindow + 1, TimedifferenceGrid1, TimedifferenceGrid2, TimedifferenceGrid3, TimedifferenceGrid4, ShortTimegrid1, sm_BGE, sm_BGN, sn_BGE, sn_BGN, sl_BGE, sl_BGN, sp_BGE, sp_BGN,  TimedifferenceGrid5, TimedifferenceGrid6,xygridsize,zgridsize);
    Resultx(i) = (X(1) - ceil(xygridsize/2))*(XYgridlimit/floor(xygridsize/2));
    Resulty(i) = (Y(1) - ceil(xygridsize/2))*(XYgridlimit/floor(xygridsize/2));
    Resultz(i) = (Z(1) - 1)*(Gridsize);
    Abstime(i) = globalstart + (i - 1) * overlapmicroseconds;
    maxoverstdxyz(i) = maxoverstd;
end

filename = sprintf('%dus_slidingwindow_%dus_overlap_%d_%d.mat', slidewindowmicroseconds, overlapmicroseconds, globalstart, globalend);
save(filename, 'Resultx', 'Resulty', 'Resultz', 'maxoverstdxyz', 'Abstime');

toc;

function [X, Y, Z,maxoverstd] = crosscorfunc(startposition, Endposition, timedifference1, timedifference2, timedifference3, timedifference4, ShortTimegrid1, sm_BGE, sm_BGN, sn_BGE, sn_BGN, sl_BGE, sl_BGN, sp_BGE, sp_BGN, timedifference5, timedifference6,xygridsize, zgridsize)
 
    %% Generate Ixy matrix
    
    Ixy = zeros(length(ShortTimegrid1), 1);
    Ixy2 = zeros(length(ShortTimegrid1), 1);
    Ixy3 = zeros(length(ShortTimegrid1), 1);
    Ixy4 = zeros(length(ShortTimegrid1), 1);
    Ixy5 = zeros(length(ShortTimegrid1), 1);
    Ixy6 = zeros(length(ShortTimegrid1), 1);
    
    ThresholdBoxlower = -0.1;
    ThresholdBoxupper = 0.1;
 
    BGE_signal = true;
    BGN_signal = false;
    BGE_Polarity = 0;
    BGN_Polarity = 0;
    buffersize = 400;
    
    % Determine which signal to use BGE or BGN based on the max(abs) value of the signal 
    
    % Choose Duke as signal polarity reference 
  
    if max(abs(sm_BGE(startposition:1:Endposition))) > max(abs(sm_BGN(startposition:1:Endposition)))
        BGE = true;
        BGN = false;
        UpperLimit = sm_BGE(startposition:1:Endposition) > ThresholdBoxupper ;
        LowerLimit = sm_BGE(startposition:1:Endposition) < ThresholdBoxlower ;
        UpperLimitsn = sn_BGE(startposition:1:Endposition+buffersize) > ThresholdBoxupper ;
        LowerLimitsn = sn_BGE(startposition:1:Endposition+buffersize) < ThresholdBoxlower ;
        UpperLimitsl = sl_BGE(startposition:1:Endposition+buffersize) > ThresholdBoxupper ;
        LowerLimitsl = sl_BGE(startposition:1:Endposition+buffersize) < ThresholdBoxlower ;
        UpperLimitsp = sp_BGE(startposition:1:Endposition+buffersize) > ThresholdBoxupper ;
        LowerLimitsp = sp_BGE(startposition:1:Endposition+buffersize) < ThresholdBoxlower ;
        
        if find(UpperLimit,1) < find(LowerLimit,1)    
            % The polarity of first pulse of BGE is positive -> All BGE must be positive 
            BGE_Polarity = 1;
            vq_sm = sm_BGE(startposition:1:Endposition); 
            
            sn_polarity = find(UpperLimitsn,1) < find(LowerLimitsn,1);
            sl_polarity = find(UpperLimitsl,1) < find(LowerLimitsl,1);
            sp_polarity = find(UpperLimitsp,1) < find(LowerLimitsp,1);
            
            if sn_polarity == 0 
                vq_sn = -1*sn_BGE(startposition:1:Endposition+buffersize);
            else
                vq_sn = 1*sn_BGE(startposition:1:Endposition+buffersize);
            end
            
            if sl_polarity == 0 
                vq_sl = -1*sl_BGE(startposition:1:Endposition+buffersize);
            else
                vq_sl = 1*sl_BGE(startposition:1:Endposition+buffersize);
            end
             
            if sp_polarity == 0 
                vq_sp = -1*sp_BGE(startposition:1:Endposition+buffersize);
            else
                vq_sp = 1*sp_BGE(startposition:1:Endposition+buffersize);
            end
            

        else
            BGE_Polarity = -1;
            % The polarity of first pulse of BGE is negative -> All BGE must be negative
            vq_sm = sm_BGE(startposition:1:Endposition); 
            sn_polarity = find(UpperLimitsn,1) < find(LowerLimitsn,1);
            sl_polarity = find(UpperLimitsl,1) < find(LowerLimitsl,1);
            sp_polarity = find(UpperLimitsp,1) < find(LowerLimitsp,1);
            
            if sn_polarity == 0 
                vq_sn = 1*sn_BGE(startposition:1:Endposition+buffersize);
            else
                vq_sn = -1*sn_BGE(startposition:1:Endposition+buffersize);
            end
            
            if sl_polarity == 0 
                vq_sl = 1*sl_BGE(startposition:1:Endposition+buffersize);
            else
                vq_sl = -1*sl_BGE(startposition:1:Endposition+buffersize);
            end
             
            if sp_polarity == 0 
                vq_sp = 1*sp_BGE(startposition:1:Endposition+buffersize);
            else
                vq_sp = -1*sp_BGE(startposition:1:Endposition+buffersize);
            end
            
 
            
        end
    else
        
        BGN = true;
        BGE = false;
        UpperLimit = sm_BGN(startposition:1:Endposition) > ThresholdBoxupper ;
        LowerLimit = sm_BGN(startposition:1:Endposition) < ThresholdBoxlower ;
        UpperLimitsn = sn_BGN(startposition:1:Endposition+buffersize) > ThresholdBoxupper ;
        LowerLimitsn = sn_BGN(startposition:1:Endposition+buffersize) < ThresholdBoxlower ;
        UpperLimitsl = sl_BGN(startposition:1:Endposition+buffersize) > ThresholdBoxupper ;
        LowerLimitsl = sl_BGN(startposition:1:Endposition+buffersize) < ThresholdBoxlower ;
        UpperLimitsp = sp_BGN(startposition:1:Endposition+buffersize) > ThresholdBoxupper ;
        LowerLimitsp = sp_BGN(startposition:1:Endposition+buffersize) < ThresholdBoxlower ;
        
        if find(UpperLimit,1) < find(LowerLimit,1) 
            BGN_Polarity = 1;
            vq_sm = sm_BGN(startposition:1:Endposition); 
            
            sn_polarity = find(UpperLimitsn,1) < find(LowerLimitsn,1);
            sl_polarity = find(UpperLimitsl,1) < find(LowerLimitsl,1);
            sp_polarity = find(UpperLimitsp,1) < find(LowerLimitsp,1);
            
            if sn_polarity == 0 
                vq_sn = -1*sn_BGN(startposition:1:Endposition+buffersize);
            else
                vq_sn = 1*sn_BGN(startposition:1:Endposition+buffersize);
            end
            
            if sl_polarity == 0 
                vq_sl = -1*sl_BGN(startposition:1:Endposition+buffersize);
            else
                vq_sl = 1*sl_BGN(startposition:1:Endposition+buffersize);
            end
             
            if sp_polarity == 0 
                vq_sp = -1*sp_BGN(startposition:1:Endposition+buffersize);
            else
                vq_sp = 1*sp_BGN(startposition:1:Endposition+buffersize);
            end

            
        else
            BGN_Polarity = -1;
            vq_sm = sm_BGN(startposition:1:Endposition);
            % The polarity of first pulse of BGE is negative -> All BGE must be negative
            sn_polarity = find(UpperLimitsn,1) < find(LowerLimitsn,1);
            sl_polarity = find(UpperLimitsl,1) < find(LowerLimitsl,1);
            sp_polarity = find(UpperLimitsp,1) < find(LowerLimitsp,1);
            
            if sn_polarity == 0 
                vq_sn = 1*sn_BGN(startposition:1:Endposition+buffersize);
            else
                vq_sn = -1*sn_BGN(startposition:1:Endposition+buffersize);
            end
            
            if sl_polarity == 0 
                vq_sl = 1*sl_BGN(startposition:1:Endposition+buffersize);
            else
                vq_sl = -1*sl_BGN(startposition:1:Endposition+buffersize);
            end
             
            if sp_polarity == 0 
                vq_sp = 1*sp_BGN(startposition:1:Endposition+buffersize);
            else
                vq_sp = -1*sp_BGN(startposition:1:Endposition+buffersize);
            end
            
            
        end
        
    end
    
 
    %% Dot product calculation based on shiftlength:
    parfor y = 1:1:length(ShortTimegrid1)
        vq_sntest = vq_sn;
        vq_sltest = vq_sl;
        vq_sptest = vq_sp;
        endpoint = 251;
        
        sn_test_shift = vq_sntest(201 + ShortTimegrid1(y):1:endpoint + ShortTimegrid1(y));
        sl_test_shift = vq_sltest(201 + ShortTimegrid1(y):1:endpoint + ShortTimegrid1(y));
        sp_test_shift = vq_sptest(201 + ShortTimegrid1(y):1:endpoint + ShortTimegrid1(y));
        
        Ixy(y) = dot(vq_sm, sn_test_shift);     %Duke Hudson
        Ixy2(y) = dot(vq_sm, sl_test_shift);    %Duke Ps2
        Ixy3(y) = dot(vq_sm, sp_test_shift);    %Duke Ps3  
        Ixy4(y) = dot(vq_sn(201:1:endpoint),sl_test_shift ); %Hudson Ps2
        Ixy5(y) = dot(vq_sn(201:1:endpoint), sp_test_shift); % Hudson Ps3
        Ixy6(y) = dot(vq_sl(201:1:endpoint), sp_test_shift); % Ps2 Ps3

        

        %TimeDiff 1 Duke -Hudson 
        %TimeDiff 2 Duke -Ps2 
        %TimeDiff 3 Duke -Ps3
        %TimeDiff 4 Hudson -Ps2
        %TimeDiff 5 Hudson -Ps3 
        %TimeDiff 6 Ps2 - Ps3 

        
    end

    %% Fill in the grids
    
    Ixyz1 =  zeros(xygridsize * xygridsize * zgridsize, 1);
    Ixyz2 =  zeros(xygridsize * xygridsize * zgridsize, 1);
    Ixyz3 =  zeros(xygridsize * xygridsize * zgridsize, 1);
    Ixyz4 =  zeros(xygridsize * xygridsize * zgridsize, 1);
    Ixyz5 =  zeros(xygridsize * xygridsize * zgridsize, 1);
    Ixyz6 =  zeros(xygridsize * xygridsize * zgridsize, 1);

    for i = 1:1:length(Ixyz1)   
        Ixyz1(i)  = Ixy(timedifference1(i) + 201);
        Ixyz2(i) = Ixy2(timedifference2(i) + 201);
        Ixyz3(i) = Ixy3(timedifference3(i) + 201);
        Ixyz4(i) = Ixy4(timedifference4(i) + 201);
        Ixyz5(i) = Ixy5(timedifference5(i) + 201);
        Ixyz6(i) = Ixy6(timedifference6(i) + 201);


    end

    %% Normalized Ixy , Sum and Plot
    
    E = Ixyz1/max(Ixyz1) + Ixyz2/max(Ixyz2) + Ixyz3/max(Ixyz3) +  Ixyz4/max(Ixyz4) + Ixyz5/max(Ixyz5) + Ixyz6/max(Ixyz6);
    
    A = reshape(E, [xygridsize xygridsize zgridsize]);

   figure(1)
   imagesc(A(:,:,20))
    
    STDImage = std(E(:));
    
    maxMatrix4 = max(E(:));
    
    maxoverstd = maxMatrix4/STDImage;
    
    [X, Y, Z] = ind2sub([xygridsize xygridsize zgridsize], find(E == maxMatrix4));

   
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
 
    Dists = R * acos(sin(lat1) .* sin(lat2) + cos(lat1) .* cos(lat2) .* cos(lon2 - lon1));
 
end

%% Time difference between random x,y coordinate to two sensors

function timedelay = calcudis(x1, y1, z1, x2, y2, z2, x, y, z)
    c = 299792.458;
    % Grid to Duke - Grid to other station
    timedelay = (sqrt((x - x1) ^ (2) + (y - y1) ^ (2) + (z - z1) ^ (2)) / c) - (sqrt((x - x2) ^ (2) + (y - y2) ^ (2) + (z - z2) ^ (2)) / c);
end

function [i ,j, k] = linearto3d(index, xgridsize, ygridsize)

i = mod(index-1,xgridsize)+1;
j = mod((index - i)/xgridsize, ygridsize)+1;
k = ((index - i)/xgridsize-j+1)/ygridsize + 1;

end
