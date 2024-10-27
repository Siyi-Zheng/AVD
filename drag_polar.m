clear
clc
close all 

kld = 15.5;
s_wet_s_ref = 5.9;

AR_temp = 8.49; % change this when the program returns the correct value
L_D_max = kld * sqrt(AR_temp / s_wet_s_ref);

e = 0.8; %oswald efficiency
e_TO_flaps= e - 0.05; %take off with flaps
e_TO_flaps_gear= e-0.1; %take off with flaps and gear
e_L_flaps = e - 0.1; %landing configuration
e_L_flaps_gear= e -0.15; %with flaps and gear

C_L_max = 2.5; % landing configuration
C_L = 1.8; % clean configuration
C_D_o = (e * pi * AR_temp) / (L_D_max * 2) ^ 2; %zero lift drag coeff
C_D_o_TO_flaps = C_D_o + 0.02; %flaps but gear up
C_D_o_TO_flaps_gear = C_D_o + 0.04; %takeoff configuration
C_D_o_L_flaps = C_D_o + 0.07; %landing configuration
C_D_o_L_flaps_gears= C_D_o + 0.09;

syms Cd
Lift_coeffcient= @(Cd, C_D_o, e) sqrt((Cd-C_D_o)/(1/(pi*AR_temp*e)));
figure
fplot(@(Cd) Lift_coeffcient(Cd, C_D_o, e), [0 0.6], 'LineWidth',1.5)

hold on
fplot(@(Cd) Lift_coeffcient(Cd, C_D_o_TO_flaps, e_TO_flaps), [0 0.6], 'LineWidth',1.5)
fplot(@(Cd) Lift_coeffcient(Cd, C_D_o_TO_flaps_gear, e_TO_flaps_gear), [0 0.6], 'LineWidth',1.5)
fplot(@(Cd) Lift_coeffcient(Cd, C_D_o_L_flaps, e_L_flaps), [0 0.6], 'LineWidth',1.5)
fplot(@(Cd) Lift_coeffcient(Cd, C_D_o_L_flaps_gears, e_L_flaps_gear), [0 0.6], 'LineWidth',1.5)
hold off
legend('Clean', 'take off with flaps','take off with flaps and gear','landing with flaps','landing with flaps and gear')
title('Drag polar')
xlabel('Cd')
ylabel('Cl')


