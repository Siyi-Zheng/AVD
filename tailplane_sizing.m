clear
clc

% typical values for airliner
V_h = 1;
V_v = 0.09;

% values from initial sizing
S_ref = 482;
b = 32.5;
mac = S_ref / (2 * b);
x_w = 30; % idk
x_cg = 25;
x_h = 70;
x_v = 65;

% sizing
S_h = V_h * S_ref * mac / (x_h - x_w);
S_v = V_v * S_ref * b / (x_v - x_w);

% pitching moment
C_m = C_m_0 + (x_cg - x_w) * C_l - eta_h * V_h * C_l_h;

% static margin
K_n = x_w - x_cg + eta_h * V_h * (1 - downwash_derivative) * a_h / a;