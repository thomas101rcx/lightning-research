clear; format short e
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

Noise = sm1.Duke_BGE(1:1:First_Pulse_index)