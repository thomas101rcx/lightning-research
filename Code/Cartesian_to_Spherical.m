clear; format short e
%% Generate rho and phi plot

load('./2017-07-08_20-58-10_Flash_50us_slidingwindow_25us_overlap_ 9.696671- 9.696721_Interp/50us_slidingwindow_25us_overlap_ 9.696671- 9.696721.mat')

[Rho, theta, phi]=CtoS(Resultx, Resulty, Resultz);



%% Convert from Cartesian coordinates back to spherical

function [rho, theta, phi] = CtoS(x, y, z)

    rho = sqrt(x.^(2) + y.^(2) + z.^(2));
    theta = acosd(z./rho);
    phi = atan2d(y,x);


end
