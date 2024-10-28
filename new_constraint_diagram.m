clc
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 = TO
% 2 = climb to cruise 1
% 3 = cruise 1 (40 000ft)
% 4 = descent to land (miss)
% 5 = climb to cruise 2 (25 000 ft)
% 6 = cruise 2
% 7 = loiter
% 8 = descent to land
% 9 = landing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MTOW AND WEIGHT FRACTIONS

%Wo = Wc + Wp + (Wf/Wo)Wo + (We/Wo)Wo
%Wo = (Wc + Wp)/(1- Wf/Wo - We/Wo)

%INITIAL WEIGHT
no_economy_pass = 416;
no_business_pass = 84;
no_crew = 12;
economy_baggage = 20;
business_baggage = 40;
crew_baggage = 10;

mass_pass = 75;

Wc = (mass_pass + crew_baggage) * no_crew;
Wp = (no_economy_pass + no_business_pass) * mass_pass + economy_baggage * no_economy_pass + business_baggage * no_business_pass;

%L/D max
kld = 15.5;
s_wet_s_ref = 5.5;

AR_temp = 8.9; % change this when the program returns the correct value
L_D_max = kld * sqrt(AR_temp / s_wet_s_ref);

%Weight fractions

sfc_cruise_mg = 13.825; %in mg/Ns
sfc_loiter_mg = 11.1;

sfc_cruise = sfc_cruise_mg * 1e-6 * 9.81; %in 1/s
sfc_loiter = sfc_loiter_mg * 1e-6 * 9.81;

%cruise 1 parameters
Range_1 = 7500 * 1852; %change from nautical miles to metres
Velocity_cruise_1 = 0.83 * 295; %speed of sound at 40,000 feet
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

figure(1)
plot(Wo1, We_Wo, "-b", LineWidth = 1.5);
hold on
plot(Wo2, We_Wo, "r", LineWidth = 1.5);
hold off
xlabel("W_o (kg)");
ylabel("W_e/W_o");
legend("W_o guess", "W_o calculated")
grid on

% function to get We from the graph (where the lines intersect)
% iterate through We_Wo and find the point where Wo1 and Wo2 are closest
close_list = We_Wo; % preallocate to the same size as We_Wo

for i = 1:length(We_Wo)
    close_list(i) = abs(Wo1(i) - Wo2(i));
end

[~, index] = min(close_list);
Wo = Wo2(index);
We_Wo_final = We_Wo(index);
We = We_Wo_final * Wo;
Wf = Wf_Wo * Wo;

figure(1)
hold on
xline(Wo, 'Color', [1 0.6 0.2], LineWidth = 1.5);
legend("W_o guess", "W_o calculated", 'Converged value')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CONSTRAINT DIAGRAM

% aircraft parameters
NE = 4; %number of engines
max_wingspan = 65; %max wingspan in m
e = 0.85; %oswald efficiency
e_TO_up = e - 0.05; % flaps but gear up
e_TO = e - 0.1; %takeoff configuration
e_L = e - 0.15; %landing configuration
C_L_max = 2.5; % landing configuration
C_L = 1.4; % clean configuration
C_D_o = (e * pi * AR_temp) / (L_D_max * 2) ^ 2; %zero lift drag coeff
C_D_o_up = C_D_o + 0.02; %flaps but gear up
C_D_o_TO = C_D_o + 0.06; %takeoff configuration
C_D_o_L = C_D_o + 0.11; %landing configuration
% cruise parameters
V_cruise_1 = 245; %cruise 1 speed in m/s
V_cruise_2 = 241.5; %cruise 2 speed in m/s (M=0.78)
rho_cruise_1 = 0.38; %cruise 1 density in kg/m^3
rho_cruise_2 = 0.54; %cruise 2 density in kg/m^3
% ground parameteres
S_L = 2900; %min TO distance / max landing distance
S_a = 305; %obstacle clearance distance (m)
% loiter parameters
rho_sl = 1.225; %density at sea level in kg/m^3
rho_ceil = 0.2372; %density at ceiling (45000ft)
rho_loiter = 1.0555; %density at loiter (5000ft)
V_loiter = 120; %velocity when loitering in m/s
% calculated values
C_L_TO = C_L_max / 1.21; % takeoff configuration
constraint_landing = ((S_L - S_a) .* C_L_max) * 9.81 / (5 * 5/3);
S_W_min = Wo * 9.81 / constraint_landing; %min reference area OF WING - seems to be around 85-90% of total area
AR = max_wingspan ^ 2 / S_W_min; %max aspect ratio
V_s = sqrt((2 * Wo * 9.81) / (rho_sl * S_W_min * C_L_max));
% landing and stall
constraint_stall = 0.5 * rho_sl * V_s * V_s * C_L_max;

% output weight Wo and AR to the console. Don't use scientific notation for weight
print_weight = floor(Wo);
fprintf("Weight: %.0f kg\n", print_weight);
fprintf("Max. wing aspect ratio: %.3g\n", AR);
disp(S_W_min)

% other values
n = [1, 1, 1, 1, 1, 1, 1, 1, 1]; %load factor for each mission stage
To = 316000 * NE; %100 % thrust for our engine in N (we need to choose one)
V_2 = 1.2 * V_s; % in takeoff configuration
V_3 = 1.25 * V_s; % used in the third phase of the climb

%thrust ratios (if it is zero it means that it is not needed for the
%calculations)
T_To = [1, (0.6531 / rho_sl) ^ 0.7, (0.244/0.2972) * (0.2972) ^ 0.7, 0, (0.833 / rho_sl) ^ 0.7, (rho_cruise_2 / rho_sl) ^ 0.7, (rho_loiter / rho_sl) ^ 0.7, 0, 0];
a_b = MWF ./ T_To; %alpha/beta = weight ratios / thrust ratios
TOP = S_L / (0.297 - 0.019 * NE);
W_S = 1:100:14000;

% CONSTRAINTS
% takeoff
constraint_TO = W_S .* (1 / (TOP * C_L_TO));

% climb rates (first climb):
% segment 1: 5% gradient, flaps to takeoff position, gear down, 3/4 engines running, speed V2
climb_rate_1a = 0.05 * V_2;
constraint_climb_1a = (4/3) * a_b(2) * (1 / (V_2)) * climb_rate_1a + ((0.5 * rho_sl * V_2 * V_2 * C_D_o_TO) ./ (MWF(2) .* W_S) + (MWF(2) * n(2) * n(2) .* W_S) ./ (0.5 * rho_sl * V_2 * V_2 * pi * AR * e_TO));
% segment 2: 3% gradient, flaps to takeoff position, gear up, 3/4 engines running, speed V2
climb_rate_1b = 0.03 * V_2;
constraint_climb_1b = (4/3) * a_b(2) * (1 / (V_2)) * climb_rate_1b + ((0.5 * rho_sl * V_2 * V_2 * C_D_o_up) ./ (MWF(2) .* W_S) + (MWF(2) * n(2) * n(2) .* W_S) ./ (0.5 * rho_sl * V_2 * V_2 * pi * AR * e_TO_up));
% segment 3: 1.7% gradient, clean configuration, 3/4 engines running, speed 1.25 * V_s (call this v3 for now)
climb_rate_1c = 0.017 * V_3;
constraint_climb_1c = (4/3) * a_b(2) * (1 / (V_3)) * climb_rate_1c + ((0.5 * rho_sl * V_3 * V_3 * C_D_o) ./ (MWF(2) .* W_S) + (MWF(2) * n(2) * n(2) .* W_S) ./ (0.5 * rho_sl * V_3 * V_3 * pi * AR * e));
% first cruise and ceilings
constraint_cruise_1 = a_b(3) .* ((0.5 * rho_cruise_1 * V_cruise_1 * V_cruise_1 * C_D_o) ./ (MWF(3) .* W_S) + (MWF(3) * n(3) * n(3) .* W_S) ./ (0.5 * rho_cruise_1 * V_cruise_1 * V_cruise_1 * pi * AR * e));
constraint_serv_ceiling = a_b(3) .* ((0.5 * rho_cruise_1 * V_cruise_1 * V_cruise_1 * ((0.85/0.83) ^ 2) * C_D_o) ./ (MWF(3) .* W_S) + (MWF(3) * n(3) * n(3) .* W_S) ./ (0.5 * rho_cruise_1 * V_cruise_1 * V_cruise_1 * ((0.85/0.83) ^ 2) * pi * AR * e));
constraint_abs_ceiling = a_b(3) .* ((0.5 * rho_ceil * V_cruise_1 * V_cruise_1 * C_D_o) ./ (MWF(3) .* W_S) + (MWF(3) * n(3) * n(3) .* W_S) ./ (0.5 * rho_ceil * V_cruise_1 * V_cruise_1 * pi * AR * e));
% climb rates (second climb):
% segment 1: 5% gradient, flaps to takeoff position, gear down, 3/4 engines running, speed V2
climb_rate2 = 0.05 * V_2;
constraint_climb_2a = (4/3) * a_b(5) * (1 / (V_2)) * climb_rate2 + ((0.5 * rho_sl * V_2 * V_2 * C_D_o_TO) ./ (MWF(5) .* W_S) + (MWF(5) * n(5) * n(5) .* W_S) ./ (0.5 * rho_sl * V_2 * V_2 * pi * AR * e_TO));
% segment 2: 3% gradient, flaps to takeoff position, gear up, 3/4 engines running, speed V2
climb_rate2 = 0.03 * V_2;
constraint_climb_2b = (4/3) * a_b(5) * (1 / (V_2)) * climb_rate2 + ((0.5 * rho_sl * V_2 * V_2 * C_D_o_up) ./ (MWF(5) .* W_S) + (MWF(5) * n(5) * n(5) .* W_S) ./ (0.5 * rho_sl * V_2 * V_2 * pi * AR * e_TO_up));
% segment 3: 1.7% gradient, clean configuration, 3/4 engines running, speed 1.25 * V_s (call this v3 for now)
climb_rate2 = 0.017 * V_3;
constraint_climb_2c = (4/3) * a_b(5) * (1 / (V_3)) * climb_rate2 + ((0.5 * rho_sl * V_3 * V_3 * C_D_o) ./ (MWF(5) .* W_S) + (MWF(5) * n(5) * n(5) .* W_S) ./ (0.5 * rho_sl * V_3 * V_3 * pi * AR * e));
% second cruise
constraint_cruise_2 = a_b(6) .* ((0.5 * rho_cruise_2 * V_cruise_2 * V_cruise_2 * C_D_o) ./ (MWF(6) .* W_S) + (MWF(6) * n(6) * n(6) .* W_S) ./ (0.5 * rho_cruise_2 * V_cruise_2 * V_cruise_2 * pi * AR * e));
% loiter
constraint_loiter = a_b(7) .* ((0.5 * rho_loiter * V_loiter * V_loiter * C_D_o) ./ (MWF(7) .* W_S) + (MWF(7) * n(7) * n(7) .* W_S) ./ (0.5 * rho_loiter * V_loiter * V_loiter * pi * AR * e));

%landing and stall and climbs stuff for plotting
T_W = [0, 1.5];
s_w = [0, 14000];
landing = [constraint_landing, constraint_landing];
stall = [constraint_stall, constraint_stall];

% get the T/W required at the landing constraint according to the takeoff constraint
T_W_landing = constraint_landing / (TOP * C_L_TO);
T_design = T_W_landing * Wo * 9.81;
T_engine = T_design / NE;
print_thrust = floor(T_engine);
fprintf("Min. thrust per engine: %.0f N\n", print_thrust);

figure(2)
% takeoff
plot(W_S, constraint_TO, color = [0.7 0 0], LineWidth = 1.5);
hold on
% climbs 1a, 1b, 1c
plot(W_S, constraint_climb_1a, color = [0 0.7 0], LineWidth = 1.5, linestyle = "-");
plot(W_S, constraint_climb_1b, color = [0 0.7 0], LineWidth = 1.5, linestyle = "--");
plot(W_S, constraint_climb_1c, color = [0 0.7 0], LineWidth = 1.5, linestyle = ":");
% cruise 1, service ceiling, absolute ceiling
plot(W_S, constraint_cruise_1, color = [0 0 0.7], LineWidth = 1.5);
plot(W_S, constraint_serv_ceiling, color = [0 0 0.7], LineWidth = 1.5, linestyle = "--");
plot(W_S, constraint_abs_ceiling, color = [0 0 0.7], LineWidth = 1.5, linestyle = ":");
% climbs 2a, 2b, 2c
plot(W_S, constraint_climb_2a, color = [0.7 0 0.7], LineWidth = 1.5, linestyle = "-");
plot(W_S, constraint_climb_2b, color = [0.7 0 0.7], LineWidth = 1.5, linestyle = "--");
plot(W_S, constraint_climb_2c, color = [0.7 0 0.7], LineWidth = 1.5, linestyle = ":");
% cruise 2, loiter
plot(W_S, constraint_cruise_2, color = [0 0.7 0.7], LineWidth = 1.5);
plot(W_S, constraint_loiter, color = [0 0.7 0.7], LineWidth = 1.5, linestyle = "--");
% landing, stall
plot(landing, s_w, color = [0.7 0.7 0], LineWidth = 1.5);
plot(stall, s_w, color = [0.7 0.7 0], LineWidth = 1.5, linestyle = "--");
scatter(constraint_landing-100, 306000*4/(Wo*9.81), "kx", "LineWidth", 1); % change this when we know our thrust and wing area
hold off
xlim([0 13000]);
ylim([0 0.6]);
grid on
xlabel("Wing Loading (kg/m^2)");
ylabel("Sea-level Thrust-to-Weight Ratio");
legend("Takeoff", "Climb 1a", "Climb 1b", "Climb 1c", "Cruise 1", "Service Ceiling", "Absolute Ceiling", "Climb 2a", "Climb 2b", "Climb 2c", "Cruise 2", "Loiter", "Landing", "Stall", "Design Point");

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
