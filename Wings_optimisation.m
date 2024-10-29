clear
clc 
close all

altitude_data= [20000, 316.0,   0.652;
                  25000, 309.6,   0.549 ;  
                  30000, 303.1,   0.459 ;     
                  350000 , 295.4,  0.38 ;
                  40000, 294.9,  0.304];

%wings optimisation
%INITIAL WEIGHT
no_economy_pass = 416;
no_business_pass = 84;
no_crew = 12;
economy_baggage = 18;
business_baggage = 40;
crew_baggage = 10;

mass_pass = 75;

Wc = (mass_pass + crew_baggage) * no_crew;
Wp = (no_economy_pass + no_business_pass) * mass_pass + economy_baggage * no_economy_pass + business_baggage * no_business_pass;

%L/D max
AR_temp_values= 7:0.1:10;
s_wet_s_ref_values= 5:0.05:6.5;
% Wo_results= zeros(size(AR_temp_values));

for i= 1:length(AR_temp_values)
    for abc= 1:length(s_wet_s_ref_values)
AR_temp = AR_temp_values(i);
kld = 15.5;
s_wet_s_ref = s_wet_s_ref_values(abc);
 % change this when the program returns the correct value
L_D_max = kld * sqrt(AR_temp / s_wet_s_ref);

%Weight fractions

sfc_cruise_mg = 13.825; %in mg/Ns
sfc_loiter_mg = 11.1;

sfc_cruise = sfc_cruise_mg * 1e-6 * 9.81; %in 1/s
sfc_loiter = sfc_loiter_mg * 1e-6 * 9.81;

cruising_altitude_index= 4;

%cruise 1 parameters
Range_1 = 7500 * 1852; %change from nautical miles to metres
Velocity_cruise_1 = 0.83 * 295.4; %speed of sound at 40,000 feet
L_D_cruise_1 = 0.866 * L_D_max;

%cruise 2 parameters
Range_2 = 370 * 1000; %change from km to m
Velocity_cruise_2 = 0.78 * 309.6; %cruise mach at 25,000 feet is 0.78
L_D_cruise_2 = 0.866 * L_D_max;

%loiter parameters
endurance = 45 * 60; %change from minutes to seconds
L_D_loiter = L_D_max;

%Let We/Wo be ef (empty fraction)
%We/Wo = ef = AWo^C

%Typical Wo of Boeing 777 = 247200 kg
%Typical Wo of A350 = 283000 kg
%Typical Wo of A340 = 276500 kg

Wo1 = 100000:25:700000;
A = 0.97;
C = -0.06;
We_Wo = (A .* (Wo1 .^ C)) * 0.95; %correction factor composite

W_fixed = Wc + Wp;

%Mission weight fractions (MWF)

MWF = MWF_Calculation(sfc_cruise, sfc_loiter, Range_1, Velocity_cruise_1, L_D_cruise_1, Range_2, Velocity_cruise_2, L_D_cruise_2, endurance, L_D_loiter);

W9_Wo = prod(MWF);
Wf_Wo = 1.015 * (1 - W9_Wo);
Wo2 = (W_fixed) ./ (1 - Wf_Wo - We_Wo);

% figure(1)
% plot(Wo1, We_Wo, "-b", LineWidth = 1.5);
% hold on
% plot(Wo2, We_Wo, "r", LineWidth = 1.5);
% hold off
% xlabel("W_o (kg)");
% ylabel("W_e/W_o");
% legend("W_o guess", "W_o calculated")
% grid on

% function to get We from the graph (where the lines intersect)
% iterate through We_Wo and find the point where Wo1 and Wo2 are closest
close_list = We_Wo; % preallocate to the same size as We_Wo

for j= 1:length(We_Wo)
    close_list(j) = abs(Wo1(j) - Wo2(j));
end

[~, index] = min(close_list);
Wo= Wo2(index);
Wo_results(i,abc)=Wo;
We_Wo_final = We_Wo(index);
We = We_Wo_final * Wo;
Wf = Wf_Wo * Wo;

% figure(1)
% hold on
% xline(Wo, 'Color', [1 0.6 0.2], LineWidth = 1.5);
% legend("W_o guess", "W_o calculated", 'Converged value')
% hold off
    end
end

[AR_temp_grid, s_wet_grid] = meshgrid(AR_temp_values, s_wet_s_ref_values);

% Plot the surface
figure;
surf(AR_temp_grid, s_wet_grid, Wo_results);
xlabel('Aspect Ratio (AR\_temp)');
ylabel('s_wet');
zlabel('wo');
title('Sensitivity of Wo to AR\_temp and s_wet_sref');
colorbar;

figure;
contourf(AR_temp_grid, s_wet_grid, Wo_results', 20); % 20 contours for more detail
xlabel('Aspect Ratio (AR\_temp)');
ylabel('Specific Fuel Consumption (sfc\_cruise)');
title('Contour Plot of Sensitivity of We to AR\_temp and sfc\_cruise');
colorbar;

% figure;
% plot(AR_temp_values, Wo_results, '-o');
% xlabel('Aspect Ratio (AR\_temp)');
% ylabel('Wo');
% title('Sensitivity of We to AR\_temp');
% grid on;

W_S_takeoff= 7900;
Sref= 482;
MTOW= 390000;
fuel_weight= 198928;
T_W_TO= 0.2898;

weight_average_cruise = (MWF(1)*MWF(2)*0.8055)*MTOW*9.81;
q_cruise= 0.5*altitude_data(4,3)*(0.83 * 294)^2;

design_cl_average= (1/q_cruise) * (weight_average_cruise/Sref)

%Function for weight fractions

function MWF = MWF_Calculation(sfc_cruise, sfc_loiter, Range_1, Velocity_cruise_1, L_D_cruise_1, Range_2, Velocity_cruise_2, L_D_cruise_2, endurance, L_D_loiter)

    MWF(1, 1) = 0.97; %fixed value
    MWF(1, 2) = 0.985; % fixed value
    MWF(1, 3) = exp(- ((Range_1 * sfc_cruise) / (Velocity_cruise_1 * L_D_cruise_1))); %cruise 1
    MWF(1, 4) = 0.99; %descend fixed
    MWF(1, 5) = 0.985; %climb fixed
    MWF(1, 6) = exp(- ((Range_2 * sfc_cruise) / (Velocity_cruise_2 * L_D_cruise_2))); %cruise 2
    MWF(1, 7) = exp(- ((endurance * sfc_loiter) / L_D_loiter)); %LOITER
    MWF(1, 8) = 0.99; % Descend fixed
    MWF(1, 9) = 0.995; %land fixed

end
