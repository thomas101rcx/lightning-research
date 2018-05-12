format short e 
[filename, pathname] = uigetfile({'*.mat';'*.*'},'File Selector');

S1 = load([pathname,filename]);
%S2 = load('./Fetched_LF/2017.07.08_20.58.14Event_on_2017-07-08_20-58_final_256_0.4_V010_speedup_denoised.txt');

Resultx = [S1.Resultx];
Resulty = [S1.Resulty];
Resultz = [S1.Resultz];
Image_MaxoverSTD = [S1.Image_maxoverstdxyz];
SNR_Duke_matrix = [S1.SNR_Duke_matrix];
SNR_Hudson_matrix = [S1.SNR_Hudson_matrix];
SNR_PS2_matrix = [S1.SNR_PS2_matrix];
SNR_PS3_matrix = [S1.SNR_PS3_matrix];
Abstime = [S1.Abstime];
Max_Index_Duke = [S1.Max_Duke_index_matrix];


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
close(h);
S2 = [VHF_time VHF_Azimuth VHF_Elevation];

%% Convert VHF Data
VHF_time = S2(:,1);
VHF_azimuth = S2(:,2);
VHF_elevation = S2(:,3);

%% Eliminate the preflash noise found in the sliding window code

SNR_Duke_matrix_Norm = normalize_var(SNR_Duke_matrix, 0, 1);
SNR_Hudson_matrix_Norm = normalize_var(SNR_Hudson_matrix, 0, 1);
SNR_PS2_matrix_Norm = normalize_var(SNR_PS2_matrix, 0, 1);
SNR_PS3_matrix_Norm = normalize_var(SNR_PS3_matrix, 0, 1);

Resultx = Resultx(SNR_Duke_matrix_Norm  > quantile(SNR_Duke_matrix_Norm, 0.6));
Resulty =  Resulty(SNR_Duke_matrix_Norm  > quantile(SNR_Duke_matrix_Norm, 0.6));
Resultz = Resultz(SNR_Duke_matrix_Norm  > quantile(SNR_Duke_matrix_Norm, 0.6));
Abstime = Abstime(SNR_Duke_matrix_Norm  > quantile(SNR_Duke_matrix_Norm, 0.6));
MaxoverSTD = Image_MaxoverSTD(SNR_Duke_matrix_Norm  > quantile(SNR_Duke_matrix_Norm, 0.6));

Resultx_filtered_NNGS = zeros(length(Resultx),1);
Resulty_filtered_NNGS = zeros(length(Resultx),1);
Resultz_filtered_NNGS = zeros(length(Resultx),1);
MaxoverSTD_filtered_NNGS = zeros(length(Resultx),1);
Abstime_filtered_NNGS = zeros(length(Resultx), 1);


A = color_jet(length(Resultx));

figure(20)
scatter3(Resultx, Resulty, Resultz,40,A(1:1:length(Resultx),:),'filled');
title('Filtered_Pre_Image_Process');
xlabel('km');
ylabel('km');
zlabel('km');
colormap jet
colorbar('Ticks', []);
caxis([9.5 10])


%% Nearest neighbor grid search with Tree implementation
% Main Tree with any points in 5 km 
% Side Tree with points outside this range

Tree_matrix = zeros(length(Resultx), 2, 'uint8');
Main_Tree_x= [];
Main_Tree_y= [];
Main_Tree_z= [];

Side_Tree_x = [];
Side_Tree_y = [];
Side_Tree_z = [];


Tree_matrix(1,1) = 1;

for i = 1:1:length(Resultx) 
    
   
  if isempty(Main_Tree_x) == 1 || isempty(Main_Tree_y) == 1 || isempty(Main_Tree_z) == 1
    Main_Tree_x = -5;
    Main_Tree_y = -10;
    Main_Tree_z = 5;
  else
    Main_Tree_x = Resultx(Tree_matrix(:,1) ~= 0);
    Main_Tree_y = Resulty(Tree_matrix(:,1) ~= 0);
    Main_Tree_z = Resultz(Tree_matrix(:,1) ~= 0);
  end
 
          
   Main_decision_Tree = 0;
   for j = 1:1:length(Main_Tree_x)

    if sqrt((Resultx(i) - Main_Tree_x(j))^(2) + (Resulty(i) -  Main_Tree_y(j))^(2) + (Resultz(i) - Main_Tree_z(j))^(2) ) < 5
       Main_decision_Tree = Main_decision_Tree + 1;
    else
       Main_decision_Tree = 0;
    end
 
   end
   
    
   if Main_decision_Tree > 0
       
       Main_Tree_x = vertcat(Main_Tree_x,Resultx(i));
       Main_Tree_y = vertcat(Main_Tree_y,Resulty(i));
       Main_Tree_z = vertcat(Main_Tree_z,Resultz(i));
       Tree_matrix(i,1) = 1;
   else
       
       Side_Tree_x = vertcat(Side_Tree_x, Resultx(i));
       Side_Tree_y = vertcat(Side_Tree_y, Resulty(i));
       Side_Tree_z = vertcat(Side_Tree_z, Resultz(i));
       Tree_matrix(i,1) = 0;
   end
   
end




for i = 1:1:length(Resultx)-1
    
    Old_location_x = Resultx(i);
    Old_location_y = Resulty(i);
    Old_location_z = Resultz(i);
    
     if sqrt( (Resultx(i+1) - Old_location_x)^(2) + (Resulty(i+1) - Old_location_y)^(2) + (Resultz(i+1) - Old_location_z)^(2) ) > 100
        %sqrt( (Resultx(i+1) - Old_location_x)^(2) + (Resulty(i+1) - Old_location_y)^(2) + (Resultz(i+1) - Old_location_z)^(2) ) >1 
        %abs(Resultx(i+1) - Old_location_x) > 1 && abs(Resulty(i+1) - Old_location_y) > 1 && abs(Resultz(i+1) - Old_location_z) > 1 
        Resultx_filtered_NNGS(i,1) = -10000;
        Resulty_filtered_NNGS(i,1) = -10000;
        Resultz_filtered_NNGS(i,1) = -10000;
        MaxoverSTD_filtered_NNGS(i,1) = -10000;
        Abstime_filtered_NNGS(i,1) = -10000;
        
    else
        
        Resultx_filtered_NNGS(i,1) = Resultx(i);
        Resulty_filtered_NNGS(i,1) = Resulty(i);
        Resultz_filtered_NNGS(i,1) = Resultz(i);
        MaxoverSTD_filtered_NNGS(i,1) = MaxoverSTD(i);
        Abstime_filtered_NNGS(i,1) = Abstime(i);
    end

   
end

Resultx_filtered_NNGS = Resultx_filtered_NNGS(Resultx_filtered_NNGS ~= -10000);
Resulty_filtered_NNGS = Resulty_filtered_NNGS(Resulty_filtered_NNGS ~= -10000);
Resultz_filtered_NNGS = Resultz_filtered_NNGS(Resultz_filtered_NNGS ~= -10000);
Abstime_filtered_NNGS = Abstime_filtered_NNGS(Abstime_filtered_NNGS ~= -10000);
MaxoverSTD_filtered_NNGS = MaxoverSTD_filtered_NNGS(MaxoverSTD_filtered_NNGS ~= -10000);

A = color_jet(length(Resultx_filtered_NNGS));

figure(2)
scatter3(Resultx_filtered_NNGS, Resulty_filtered_NNGS, Resultz_filtered_NNGS,40,A(1:1:length(Resultz_filtered_NNGS),:),'filled');
title('Filtered_NNGS');
xlabel('km');
ylabel('km');
zlabel('km');
colormap jet
colorbar('Ticks', []);
caxis([9.5 10])


Resultx_filtered = zeros(length(Resultx_filtered_NNGS),1);
Resulty_filtered = zeros(length(Resulty_filtered_NNGS),1);
Resultz_filtered = zeros(length(Resultz_filtered_NNGS),1);
MaxoverSTD_filtered = zeros(length(MaxoverSTD_filtered_NNGS),1);
Abstime_filtered = zeros(length(Abstime_filtered_NNGS), 1);

%% Filter the rest of the window using MaxoverSTD

for i = 1:1:length(Resultx_filtered_NNGS) 
    
    if MaxoverSTD_filtered_NNGS(i) > quantile(MaxoverSTD_filtered_NNGS, 0.6)
        Resultx_filtered(i, 1) = Resultx_filtered_NNGS(i);
        Resulty_filtered(i, 1) = Resulty_filtered_NNGS(i);
        Resultz_filtered(i, 1) = Resultz_filtered_NNGS(i);
        Abstime_filtered(i, 1) = Abstime_filtered_NNGS(i);
    else
        Resultx_filtered(i, 1) = -100000;
        Resulty_filtered(i, 1) = -100000;
        Resultz_filtered(i, 1) = -100000;
        Abstime_filtered(i, 1) = -100000;
    end
end


%% Eliminate the windows that have MaxoverSTD less than the average of all MaxoverSTD

Resultx_filtered = Resultx_filtered(Resultx_filtered ~= -100000);
Resulty_filtered = Resulty_filtered(Resulty_filtered ~= -100000);
Resultz_filtered = Resultz_filtered(Resultz_filtered ~= -100000);
Abstime_filtered = Abstime_filtered(Abstime_filtered ~= -100000);


%% Convert Cart to Sph 

[azimuth, elevation, r]=cart2sph(Resultx_filtered, Resulty_filtered, Resultz_filtered);


%% Convert Radians to degree
LF_azimuthdegree = (azimuth*360)/(2*pi);

LF_elevationdegree = (elevation*360)/(2*pi);


%% Convert azimuth to geographic azimuth

LF_azimuthdegree_geo = -1*(LF_azimuthdegree) + 90;

LF = [Abstime_filtered LF_azimuthdegree_geo LF_elevationdegree r];
VHF = [VHF_time VHF_Azimuth VHF_Elevation];
save LF.mat LF
save VHF.mat VHF


% figure(1)
% pol = polarscatter((azimuth*-1) + pi/2, elevation);
% ax = gca;
% ax.RDir = 'reverse';
% ax.ThetaZeroLocation = 'top';
% ax.ThetaDir = 'clockwise';
% ax.ThetaLim = [0 360];
% 

A = color_jet(length(Resultx_filtered));

figure(13)
scatter3(Resultx_filtered, Resulty_filtered, Resultz_filtered,40,A(1:1:length(Resultz_filtered),:),'filled');
title('Filtered_NNGS_Image_maxoverSTD');
xlabel('km');
ylabel('km');
zlabel('km');
colormap jet
colorbar('Ticks', []);
caxis([9.5 10])


figure(3)
scatter(VHF_time, VHF_azimuth);
title('VHF azimuth time');

figure(4)
scatter(VHF_time, VHF_elevation);
title('VHF elevation time');

figure(5)
scatter(Abstime_filtered, LF_azimuthdegree_geo);
title('LF azimuth time');


figure(6)
scatter(Abstime_filtered, LF_elevationdegree);
title('LF elevation time');

figure(7)
scatter(Abstime_filtered, Resultx_filtered);
title('LF X axis vs time');

figure(8)
scatter(Abstime_filtered, Resulty_filtered);
title('LF Y axis vs time');

figure(9)
scatter(Abstime_filtered, Resultz_filtered);
title('LF Z height');


figure(10)
scatter(Resultx_filtered, Resulty_filtered,25, A(1:1:length(Resultx_filtered),:),'filled')
title('LF plane view ')


figure(11)
scatter(LF_elevationdegree, Resultz_filtered)
title('LF elevation vs height');


%% Interpolate the elevation vs height for LF
%% Use this information to match VHF data and find the height for VHF -> Get rho from this information


% A = [LF_elevationdegree Resultz_filtered];
% B = unique(A,'rows','stable');
% y_VHF = interp1(B(:,1), B(:,2),VHF_Elevation(24984:1:25594), 'linear');
% 
% figure(12)
% scatter(VHF_Elevation(24984:1:25594), y_VHF);
% title('VHF elevation vs height after interpolation');
% 
% 
% VHF_rho = y_VHF./cosd(VHF_Elevation(24984:1:25594));
% 
% figure(13)
% scatter(VHF_time(24984:1:25594), VHF_rho);
% 
% % Convert from spherical back to Cart
% 
% VHF_azimuth = ((VHF_azimuth-90)*-1)*((2*pi)/360);
% 
% VHF_elevation = (VHF_elevation)*((2*pi)/360);
% 
% [VHF_x,VHF_y,VHF_z] = sph2cart(VHF_azimuth(24984:1:25594), VHF_elevation(24984:1:25594), VHF_rho);
% 
% 
% A = color_jet(length(VHF_x));
% 
% 
% figure(20)
% scatter3(VHF_x, VHF_y, VHF_z,40,A(1:1:length(VHF_x),:),'filled');
% title('Filtered');
% xlabel('km');
% ylabel('km');
% zlabel('km');
% colormap jet
% colorbar('Ticks', []);
% caxis([9.5 10])




%% Color bar matrix generation

function color_jet_arr = color_jet(size_array)

        step = ceil(size_array/3);
        color_jet_arr = zeros(size_array, 3);
        color_jet_arr(1:1:step, 1) = zeros(step, 1);
        color_jet_arr(1:1:step, 2) = 1/step:1/step:1;
        color_jet_arr(1:1:step, 3) = zeros(step, 1)+1;
        color_jet_arr(step+1:1:step*2, 1) = 1/step:1/step:1;
        color_jet_arr(step+1:1:step*2, 2) = zeros(step, 1)+1;
        color_jet_arr(step+1:1:step*2, 3) = 1:-1/step:1/step;
        color_jet_arr(2*step+1:1:step*3, 1) = zeros(step, 1)+1;
        color_jet_arr(2*step+1:1:step*3, 2) = 1:-1/step:1/step;
        color_jet_arr(2*step+1:1:step*3, 3) = zeros(step, 1);
        
end

function normalized = normalize_var(array, x, y)

     % Normalize to [0, 1]:
     m = min(array);
     range = max(array) - m;
     array = (array - m) / range;

     % Then scale to [x,y]:
     range2 = y - x;
     normalized = (array*range2) + x;
end
