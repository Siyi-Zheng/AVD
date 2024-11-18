clear
clc

% typical values for airliner
V_h = 0.535;
V_v = 0.062;

% values from initial sizing
S_ref = 482;
b = 32.5;
mac = S_ref / (2 * b);
x_w = 38;
x_cg = 37;
x_h = 75;
x_v = 73;
C_m_0 = -0.1433; % check this but its what xfoil gives
eta_h = 0.9; % typical value
downwash_derivative = 0.1; % typical value
a = 7.74;
a_h = 6; % this changes when we choose an aerofoil

% sizing
S_h = V_h * S_ref * mac / (x_h - x_w);
S_v = V_v * S_ref * b / (x_v - x_w);

% pitching moment
% C_m = C_m_0 + (x_cg - x_w) * C_l - eta_h * V_h * C_l_h;

% static margin
K_n = (x_w - x_cg + eta_h * V_h * (1 - downwash_derivative) * a_h / a) / mac;