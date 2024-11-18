%% Specific Excess Power Plot
clear;
clc;
close all;

% Constants
g = 9.81;                    % Acceleration due to gravity [m/s^2]
R = 287;                     % Gas constant for air [J/(kg*K)]
gamma = 1.4;                 % Ratio of specific heats for air

% Aircraft parameters
Ne = 4;                      % Number of engines
BPR = 7.4;                   % Bypass Ratio
AR = 8.7;                    % Aspect Ratio
e = 0.8;                     % Oswald Efficiency
W = 350000 * g;              % Weight [N]
S = 130;                     % Wing reference area in square meters
CD0 = 0.0162;                % Zero-lift drag coefficient
CDi = 0.012;
K = 1 / (pi * AR * e);       % Induced drag factor
T_sl = Ne * 258.34 * 1000;   % Sea-level static thrust [N]

% Define the grid for Mach number and altitude
mach = linspace(0, 1.2, 50);         % Mach number range
altitude = linspace(0, 10000, 50);   % Altitude range
[MACH, ALT] = meshgrid(mach, altitude);

% % Atmospheric model
% sigma_h = exp(-ALT / 7000);                % Simplified density ratio
% Temp = 288.15 - 0.0065 * ALT;              % Temperature [K]
% rho = 1.225 * sigma_h;                     % Air density [kg/m^3]
% a = sqrt(gamma * R * Temp);                % Speed of sound [m/s]
% V = MACH .* a;                             % True airspeed [m/s]
% q = 0.5 * rho .* V.^2;                     % Dynamic pressure

% used Internation Standard Atmosphere model from matlab
[T, a, P, rho] = atmosisa(ALT);
V = MACH .* a; 
q = 0.5 * rho .* V.^2;

% Lift coefficient for level flight
CL = (W ./ (q * S));                       % Required lift coefficient

% Drag Model
[D, CD] = DragModel(CD0, K, CL, q, S, MACH);

% Thrust Lapse Model
tau = ThrustLapse(ALT, MACH, BPR); 
T = T_sl * tau;                            % Available thrust [N]

% Specific Excess Power Calculation
Ps = (T - D) .* V / W;                     % Specific excess power

% Plotting Ps Contours
figure;
[c, h] = contour(MACH, ALT, Ps, 0:5:200, 'ShowText', 'on');
clabel(c, h, 'LabelSpacing', 300);
hold on;

% Add stall speed line
CL_max = 1.4;                                % Max lift coefficient
stall_speed = sqrt((2 * W) ./ (rho .* S .* CL_max));
stall_mach = stall_speed ./ a;
L1 = plot(stall_mach, ALT, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Stall Speed');

% Add cruise Mach line
cruise_mach = 0.83;                           % Cruise Mach number
L2 = plot(cruise_mach * ones(size(ALT)), ALT, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Cruise Mach');

% Add max Mach number line
max_mach = 0.9;                              % Max Mach number
L3 = plot(max_mach * ones(size(ALT)), ALT, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Max Mach Requirement');

% Labels and legend
xlabel('Mach Number',FontSize=14);
ylabel('Altitude (m)',FontSize=14);
title('Specific Excess Power (Ps) Contour',FontSize=16);
lgd = legend([h,L1(1),L2(1),L3(1)],'Ps Contours', 'Stall Speed', 'Cruise Mach', 'Max Mach Requirement', 'Location', 'best');
lgd.FontSize = 14;
grid on;

%% Drag Model Function
function [D, CD] = DragModel(CD0, K, CL, q, S, MACH)
    % Induced drag
    CDi = K * CL.^2;

    % Wave drag starts beyond Mach 0.75
    CDw = zeros(size(MACH));
    CDw(MACH > 0.75) = 0.002 * (MACH(MACH > 0.75) - 0.75).^2;

    % Total drag coefficient
    CD = CD0 + CDi + CDw;

    % Total drag force
    D = q .* CD .* S;
end

%% Thrust Lapse Function
function tau = ThrustLapse(h, M, BPR)
    % Thrust lapse for high BPR engines
    K1 = 0.89;
    K2 = -0.014;
    K3 = -0.3;
    K4 = 0.005;

    % Atmospheric density ratio
    sigma_h = exp(-h / 7000);

    % Altitude and Mach components
    B1 = zeros(size(h));
    B1(h <= 11000) = sigma_h(h <= 11000).^0.7;
    B1(h > 11000) = 1.439 * sigma_h(h > 11000);

    B2 = K1 + K2 * BPR + (K3 + K4 * BPR) .* M;

    % Thrust lapse ratio
    tau = B1 .* B2;
end
