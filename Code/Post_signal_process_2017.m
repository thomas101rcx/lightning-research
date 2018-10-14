format short e 
[filename, pathname] = uigetfile({'*.mat';'*.*'},'File Selector');

S1 = load([pathname,filename]);
%S2 = load('./Fetched_LF/2017.07.08_20.58.14Event_on_2017-07-08_20-58_final_256_0.4_V010_speedup_denoised.txt');

Resultx = [S1.Resultx];
Resulty = [S1.Resulty];
Resultz = [S1.Resultz];
Image_SNR = [S1.Image_maxoverstdxyz];
Image_max = [S1.Image_max];
Image_std = [S1.Image_std];
SNR_Duke = [S1.SNR_Duke_matrix];
SNR_Hudson = [S1.SNR_Hudson_matrix];
SNR_PS2 = [S1.SNR_PS2_matrix];
SNR_PS3 = [S1.SNR_PS3_matrix];
Abstime = [S1.Abstime];
Max_Index_Duke = [S1.Max_Duke_index_matrix];
Max_Index_Hudson = [S1.Max_Hudson_index_matrix];
Max_Index_PS2 = [S1.Max_PS2_index_matrix];
Max_Index_PS3 = [S1.Max_PS3_index_matrix];
XCorr_Duke_Hudson_max = [S1.XCorr_Duke_Hudson_matrix];
XCorr_Duke_PS2_max = [S1.XCorr_Duke_PS2_matrix];
XCorr_Duke_PS3_max = [S1.XCorr_Duke_PS3_matrix];
XCorr_Hudson_PS2_max = [S1.XCorr_Hudson_PS2_matrix];
XCorr_Hudson_PS3_max = [S1.XCorr_Hudson_PS3_matrix];
XCorr_PS2_PS3_max = [S1.XCorr_PS2_PS3_matrix];
Xcorrsum = XCorr_Duke_Hudson_max+  XCorr_Duke_PS2_max + XCorr_Duke_PS3_max + XCorr_Hudson_PS2_max + XCorr_Hudson_PS3_max + XCorr_PS2_PS3_max;

%% Filter Interpolated Signal using SNR 
% Rules : 1. Signal SNR -> Duke >  6.54 , Hudson > 4.07, PS2 > 2.36, PS3 > 2.57 , Base on quantile 50% 

SNR_quantile = 0.5;
Index_upper_Edge = 850;
Index_lower_Edge = 150;
XCorr_Coeff_quantile = 0.5;
Image_SNR_quantile = 0.5;


Image_max_matrix = [Xcorrsum, Image_max];

Error = abs(Image_max_matrix(:,1) - Image_max_matrix(:,2))./ Image_max_matrix(:,2);

Error_bar = 0.1;
Image_SNR_bar = 11 ;


% Resultx = Resultx (Error < Error_bar & Image_SNR > Image_SNR_bar);
% Resulty = Resulty (Error < Error_bar & Image_SNR > Image_SNR_bar);
% Resultz = Resultz (Error < Error_bar & Image_SNR > Image_SNR_bar);

% Resultx = Resultx (Xcorrsum > quantile(Xcorrsum, XCorr_Coeff_quantile) );
% Resulty = Resulty (Xcorrsum > quantile(Xcorrsum, XCorr_Coeff_quantile) );
% Resultz = Resultz (Xcorrsum > quantile(Xcorrsum, XCorr_Coeff_quantile) );

Resultx = Resultx( (SNR_Duke > quantile(SNR_Duke, SNR_quantile)) & (Max_Index_Duke > Index_lower_Edge & Max_Index_Duke < Index_upper_Edge) & (Xcorrsum > quantile(Xcorrsum, XCorr_Coeff_quantile)) & (Image_SNR > quantile(Image_SNR, Image_SNR_quantile)) );
Resulty = Resulty( (SNR_Duke > quantile(SNR_Duke, SNR_quantile)) & (Max_Index_Duke > Index_lower_Edge & Max_Index_Duke < Index_upper_Edge) & (Xcorrsum > quantile(Xcorrsum, XCorr_Coeff_quantile)) & (Image_SNR > quantile(Image_SNR, Image_SNR_quantile)) );
Resultz = Resultz( (SNR_Duke > quantile(SNR_Duke, SNR_quantile)) & (Max_Index_Duke > Index_lower_Edge & Max_Index_Duke < Index_upper_Edge) & (Xcorrsum > quantile(Xcorrsum, XCorr_Coeff_quantile)) & (Image_SNR > quantile(Image_SNR, Image_SNR_quantile)) );


% Max_Index_Duke = Max_Index_Duke( SNR_Duke > quantile(SNR_Duke, SNR_quantile));
% Max_Index_Hudson = Max_Index_Hudson( SNR_Duke > quantile(SNR_Duke, SNR_quantile));
% Max_Index_PS2 = Max_Index_PS2( SNR_Duke > quantile(SNR_Duke, SNR_quantile));
% Max_Index_PS3 = Max_Index_PS3( SNR_Duke > quantile(SNR_Duke, SNR_quantile));
% Image_SNR = Image_SNR(SNR_Duke > quantile(SNR_Duke, SNR_quantile));
% XCorr_Duke_Hudson_max = XCorr_Duke_Hudson_max(SNR_Duke > quantile(SNR_Duke, SNR_quantile));

figure(1) 
scatter3(Resultx,Resulty,Resultz);

%% Filter out edge signals ( < 150 || > 850)

% Resultx  = Resultx((Max_Index_Duke > Index_lower_Edge & Max_Index_Duke < Index_upper_Edge));
% Resulty  = Resulty((Max_Index_Duke > Index_lower_Edge & Max_Index_Duke < Index_upper_Edge));
% Resultz  = Resultz((Max_Index_Duke > Index_lower_Edge & Max_Index_Duke < Index_upper_Edge));
% Image_SNR = Image_SNR((Max_Index_Duke > Index_lower_Edge & Max_Index_Duke < Index_upper_Edge));
% XCorr_Duke_Hudson_max = XCorr_Duke_Hudson_max((Max_Index_Duke > Index_lower_Edge & Max_Index_Duke < Index_upper_Edge));
% 
% 
% %figure(2) 
% %scatter3(Resultx,Resulty,Resultz);
% 

%%  Filter out XCorr_Coefficient 
% 
% 
% 
% Resultx = Resultx ( XCorr_Duke_Hudson_max > quantile(XCorr_Duke_Hudson_max, XCorr_Coeff_quantile));
% Resulty = Resulty ( XCorr_Duke_Hudson_max > quantile(XCorr_Duke_Hudson_max, XCorr_Coeff_quantile));
% Resultz = Resultz ( XCorr_Duke_Hudson_max > quantile(XCorr_Duke_Hudson_max, XCorr_Coeff_quantile));
% Image_SNR = Image_SNR(XCorr_Duke_Hudson_max > quantile(XCorr_Duke_Hudson_max, XCorr_Coeff_quantile));
% 
% 
% %figure(3) 
% %scatter3(Resultx,Resulty,Resultz);
% 
%% Filter out Image_SNR
% 
% 
% 
% Resultx =  Resultx(Image_SNR > quantile(Image_SNR, Image_SNR_quantile));
% Resulty =  Resulty(Image_SNR > quantile(Image_SNR, Image_SNR_quantile));
% Resultz =  Resultz(Image_SNR > quantile(Image_SNR, Image_SNR_quantile));
% 
% 
% figure(4) 
% scatter(Resultx,Resulty);
% xlim([-15 15]);







