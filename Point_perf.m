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
AR = 8.77;                    % Aspect Ratio
e = 0.85;                     % Oswald Efficiency
W = 390000 * g;              % Weight [N]
S = 482;                     % Wing reference area in square meters
CD0 = 0.0161;                % Zero-lift drag coefficient
K = 1 / (pi * AR * e);       % Induced drag factor
T_sl = Ne * 295.8 * 1000;   % Sea-level static thrust [N]

% Define the grid for Mach number and altitude
mach = linspace(0, 1.2, 50);         % Mach number range
altitude = linspace(0, 15000, 50);   % Altitude range
[MACH, ALT] = meshgrid(mach, altitude);

% used Internation Standard Atmosphere model from matlab
[T, a, P, rho] = atmosisa(ALT);
V = MACH .* a; 
q = 0.5 * rho .* V.^2;

% Lift coefficient for level flight
CL = (W ./ (q * S));                       % Required lift coefficient

% Drag Model
[D, CD] = DragModel(CD0, K, CL, q, S, MACH);

% Thrust Lapse Model
[tau,B1,B2] = ThrustLapse(ALT, MACH, BPR); 
% tau = ThrustLapseSimple(ALT,0.83);
T = T_sl * tau;                            % Available thrust [N]

% Specific Excess Power Calculation
Ps = (T - D) .* V / W;                     % Specific excess power

% Plotting Ps Contours
figure;
[c, h] = contour(MACH, ALT * 3.2808, Ps, 0:2:200, 'ShowText', 'on');
clabel(c, h, 'LabelSpacing', 300);
hold on;

% Add stall speed line
CL_max = 1.4;                                % Max lift coefficient
stall_speed = sqrt((2 * W) ./ (rho .* S .* CL_max));
stall_mach = stall_speed ./ a;
L1 = plot(stall_mach, ALT * 3.2808, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Stall Speed');

% Add cruise Mach line
cruise_mach = 0.83;                           % Cruise Mach number
L2 = plot(cruise_mach * ones(size(ALT)), ALT * 3.2808, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Cruise Mach');

% Add max Mach number line
max_mach = 0.85;                              % Max Mach number
L3 = plot(max_mach * ones(size(ALT)), ALT * 3.2808, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Max Mach Requirement');

% Labels and legend
xlabel('Mach Number',FontSize=14);
ylabel('Altitude (ft)',FontSize=14);

lgd = legend([h,L1(1),L2(1),L3(1)],'P_s Contours', 'Stall Speed', 'Cruise Mach', 'Max Mach Requirement', 'Location', 'best');
lgd.FontSize = 14;
grid on;
hold off




%% Drag Model Function
function [D, CD] = DragModel(CD0, K, CL, q, S, MACH)
    % Induced drag
    CDi = K * CL.^2;

    % Wave drag starts beyond Mach 0.75
    CDw = zeros(size(MACH));

    % wave drag params
    l = 77.6;
    Ko = 1.5; % for an airliner
    Kf = 1.5;
    A = 8.77;
    S = 482;
    Bstar = 6.34 * 0.785;
    H = 6.34;
    Ap = 4;
    d = 1.13 * (Bstar * H - Ap) ^ 0.5;
    % equation
    CDw = Ko * (l/d) ^ (-2) * (9.4 * Kf / ((l/d)^2 * (S/l^2)) +...
        1.2 * (Kf*S/(A*l^2)) ^ 0.5 * (0.14/0.05)) * (1 + 0.0034 * (3-MACH).^3.5);
    
    % Total drag coefficient
    CD = CD0 + CDi + CDw;

    % Total drag force
    D = q .* CD .* S;
end

%% Thrust Lapse Function
function [tau,B1,B2] = ThrustLapse(h, M, BPR)
    % Thrust lapse for high BPR engines
    K1 = 0.89;
    K2 = -0.014;
    K3 = -0.3;
    K4 = 0.005;

    % Atmospheric density ratio
    % sigma_h = exp(-h / 7000);
    [~, ~, ~, rho] = atmosisa(h);
    sigma_h = rho/1.225;

    % Altitude and Mach components
    B1 = zeros(size(h));
    B1(h <= 11000) = sigma_h(h <= 11000).^0.7;
    B1(h > 11000) = 1.439 * sigma_h(h > 11000);

    B2 = K1 + K2 * BPR + (K3 + K4 * BPR) .* M;

    % Thrust lapse ratio
    tau = B1 .* B2;
end

    
    
