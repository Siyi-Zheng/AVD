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

Wc = 1020; %kg
Wp = 60000; %kg

%Let We/Wo be ef (empty fraction)
%We/Wo = ef = AWo^C

%Typical Wo of Boeing 777 = 247200 kg
%Typical Wo of A350 = 283000 kg
%Typical Wo of A340 = 276500 kg

Wo1 = 100000:25:700000;
A = 0.97;
C = -0.06;
We_Wo = A .* (Wo1 .^ C);

W_fixed = Wc + Wp;
% W_fixed = 65000;

%Mission weight fractions (MWF)
MWF = [0.97
       0.985
       0.637855042
       0.99
       0.985
       0.987929129
       0.985645559
       0.99
       0.995]';
W9_Wo = prod(MWF);
Wf_Wo = 1.015 * (1 - W9_Wo);
Wo2 = (W_fixed) ./ (1 - Wf_Wo - We_Wo);

figure(1)
plot(Wo1, We_Wo, "-b", LineWidth = 1.5);
hold on
plot(Wo2, We_Wo, "r", LineWidth = 1.5);
hold off
xlabel("W_o");
ylabel("W_e/W_o");
title("Convergence of MTOW");
legend("W_o guess", "W_o calculated")
grid on

%values taken from the plot + calculations
Wo = 410000;

% function to get We from the graph (where the lines intersect)
% iterate through We_Wo and find the point where Wo1 and Wo2 are closest
close_list = We_Wo; % preallocate to the same size as We_Wo
for i = 1:length(We_Wo)
    close_list(i) = abs(Wo1(i) - Wo2(i));
end
[~, index] = min(close_list);
We = We_Wo(index) * Wo;
Wf = Wf_Wo * Wo;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CONSTRAINT DIAGRAM

% aircraft parameters
NE = 4; %number of engines
AR = 8; %Aspect ratio
max_wingspan = 65; %max wingspan in m
e = 0.8; %oswald efficiency
C_L_max = 2.5; %max lift coeff 2:0.1:3
C_L = 1.8; %CL without flaps
L_D_max = 19.14375; %max lift to drag ratio
C_D_o = 0.017358852; %zero lift drag coeff
S_ref = 443.5; %reference area in m^2
% cruise parameters
V_cruise_1 = 245; %cruise 1 speed in m/s
V_cruise_2 = 241.5; %cruise 2 speed in m/s (M=0.78)
rho_cruise_1 = 0.3016; %cruise 1 density in kg/m^3
rho_cruise_2 = 0.54; %cruise 2 density in kg/m^3
% ground parameteres
rho_sl = 1.225; %density at sea level in kg/m^3
rho_ceil = 0.2372; %density at ceiling (45000ft)
n = [1, 1, 1, 1, 1, 1, 1, 1, 1]; %load factor for each mission stage
rho_loiter = 1.0555; %density at loiter (5000ft)
V_loiter = 120; %velocity when loitering in m/s
C_L_TO = C_L_max / 1.21; %lift coeff at TO
S_L = 2900; %min TO distance / max landing distance
S_a = 305; %obstacle clearance distance (m)
To = 316000 * 4; %100 % thrust for our engine in N (we need to choose one)
v_x = 147.8;
% calculated values
V_s = sqrt((2 * Wo * 9.81) / (rho_sl * S_ref * C_L_max));

%thrust ratios (if it is zero it means that it is not needed for the
%calculations)
T_To = [1, (0.6531 / rho_sl) ^ 0.7, (0.244/0.2972) * (0.2972) ^ 0.7, 0, (0.833 / rho_sl) ^ 0.7, (rho_cruise_2 / rho_sl) ^ 0.7, (rho_loiter / rho_sl) ^ 0.7, 0, 0];
a_b = MWF ./ T_To; %alpha/beta = weight ratios / thrust ratios
TOP = S_L / (0.297 - 0.019 * NE);
W_S = 1:100:11000;
V_s_changing = sqrt((2 * 9.81 * W_S) / (rho_sl * C_L_max));

% TODO: add the cruise equations to climbing as drag still needs to be overcome
% TODO: replace the climb rates with regulations according to the three-phase flight

constraint_TO = W_S .* (1 / (TOP * C_L_TO));
constraint_climb_1 = a_b(2) * (1 ./ (1.1 * V_s_changing)) * climb_rate1;
constraint_climb_oei = 2 .* a_b(2) .* ((0.5 * rho_cruise_1 * V_cruise_1 * V_cruise_1 * C_D_o) ./ (MWF(2) .* W_S) + (MWF(2) * n(2) * n(2) .* W_S) ./ (0.5 * rho_cruise_1 * V_cruise_1 * V_cruise_1 * pi * AR * e));
constraint_cruise_1 = a_b(3) .* ((0.5 * rho_cruise_1 * V_cruise_1 * V_cruise_1 * C_D_o) ./ (MWF(3) .* W_S) + (MWF(3) * n(3) * n(3) .* W_S) ./ (0.5 * rho_cruise_1 * V_cruise_1 * V_cruise_1 * pi * AR * e));
constraint_climb_2 = a_b(5) * (1 ./ (1.1 * V_s_changing)) * climb_rate2;
constraint_cruise_2 = a_b(6) .* ((0.5 * rho_cruise_2 * V_cruise_2 * V_cruise_2 * C_D_o) ./ (MWF(6) .* W_S) + (MWF(6) * n(6) * n(6) .* W_S) ./ (0.5 * rho_cruise_2 * V_cruise_2 * V_cruise_2 * pi * AR * e));
constraint_loiter = a_b(7) .* ((0.5 * rho_loiter * V_loiter * V_loiter * C_D_o) ./ (MWF(7) .* W_S) + (MWF(7) * n(7) * n(7) .* W_S) ./ (0.5 * rho_loiter * V_loiter * V_loiter * pi * AR * e));
constraint_landing = (((S_L * 3/5) - S_a) .* C_L_max) * 9.81 / 5;
constraint_serv_ceiling = a_b(3) .* ((0.5 * rho_cruise_1 * V_cruise_1 * V_cruise_1 * ((0.85/0.83) ^ 2) * C_D_o) ./ (MWF(3) .* W_S) + (MWF(3) * n(3) * n(3) .* W_S) ./ (0.5 * rho_cruise_1 * V_cruise_1 * V_cruise_1 * ((0.85/0.83) ^ 2) * pi * AR * e));
constraint_abs_ceiling = a_b(3) .* ((0.5 * rho_ceil * V_cruise_1 * V_cruise_1 * C_D_o) ./ (MWF(3) .* W_S) + (MWF(3) * n(3) * n(3) .* W_S) ./ (0.5 * rho_ceil * V_cruise_1 * V_cruise_1 * pi * AR * e));
constraint_stall = 0.5 * rho_sl * V_s * V_s * C_L_max;
constraint_AR = Wo * AR / (max_wingspan ^ 2) * 9.81;
constraint_TWR = To / (Wo * 9.81);

%landing and stall and climbs stuff for plotting
T_W = [0, 1.5];
s_w = [0, 11000];
climb1 = [constraint_climb_1, constraint_climb_1];
climb2 = [constraint_climb_2, constraint_climb_2];
landing = [constraint_landing, constraint_landing];
stall = [constraint_stall, constraint_stall];
aspect_ratio = [constraint_AR, constraint_AR];
thrust_weight_ratio = constraint_TWR * W_S ./ W_S;

figure(2)
plot(W_S, constraint_TO, "r", LineWidth = 1.5);
hold on
plot(s_w, climb1, LineWidth = 1.5, linestyle="--");
plot(W_S, constraint_climb_oei, LineWidth = 1.5);
plot(W_S, constraint_cruise_1, LineWidth = 1.5);
plot(s_w, climb2, LineWidth = 1.5);
plot(W_S, constraint_cruise_2, LineWidth = 1.5);
plot(W_S, constraint_loiter, LineWidth = 1.5, linestyle=":");
plot(landing, T_W, LineWidth = 1.5, linestyle="--");
plot(W_S, constraint_serv_ceiling, LineWidth = 1.5, linestyle="-.");
plot(W_S, constraint_abs_ceiling, LineWidth = 1.5, linestyle="--");
plot(stall, T_W, LineWidth = 1.5, linestyle="--");
plot(W_S, thrust_weight_ratio, LineWidth = 1.5, color="black", linestyle=":");
%plot(aspect_ratio, T_W, LineWidth = 1.5, color="black", linestyle="--");
hold off
title("Constraint Diagram");
xlim([0 10000]);
ylim([0 1]);
grid on
xlabel("W/S");
ylabel("(T/W)_o");
legend("Takeoff", "Climb 1", "Climb with one engine out", "Cruise 1",...
     "Climb 2", "Cruise 2", "Loiter", "Landing", "Service ceiling", ...
     "Absolute ceiling", "Stall", "Aircraft min T/W", "Aspect ratio = 8");
