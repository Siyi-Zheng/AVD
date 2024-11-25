clear
clc

% typical values for airliner
V_h = 0.535;
V_v = 0.062;

% values from initial sizing
S_ref = 482;
b = 32.5;
mac = S_ref / (2 * b);
x_w = 41.5;
x_cg = 39.5;
x_h = 75;
x_v = 73;
C_m_0 = -0.1433; % check this but its what xfoil gives
eta_h = 0.9; % typical value
downwash_derivative = 0.1; % typical value
a = 5.581;
a_h = 4.39; % this changes when we choose an aerofoil

% sizing
S_h = V_h * S_ref * mac / (x_h - x_w);
S_v = V_v * S_ref * b / (x_v - x_w);