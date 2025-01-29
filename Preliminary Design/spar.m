clear all
clc

% spar

a = 2; % Web panel spacing, m
E = 70; % Young's Modulus, GPa
sigma_y = 324; % Yield Stress, MPa
c = 8.06; % Wing box width, m
b2 = 1.51; % Wing box height, m
Ks = 10.4; % Read from graph
V = 3311650; % Shear Load, N
T = 303353; % Torque Load, Nm
q0 = T / (2 * c * b2 * 1000); % Torque shear flow, N/mm
q2 = V / (2 * b2 * 1000); % Load shear flow, N/mm
q_FS = q2 + q0; % Front spar shear flow, N/mm
q_RS = q2 - q0; % Rear spar shear flow, N/mm
t_FS = (q_FS * 1000 * b2 / (Ks * E * 1e9)) ^ (1/3) * 1000; % Front web thickness, mm
t_RS = (q_RS * 1000 * b2 / (Ks * E * 1e9)) ^ (1/3) * 1000; % Rear web thickness, mm

t2 = 9.6; % skin thickness, mm (FROM SKIN)
tau_0 = q0 / t2; % shear stress, N/mm
b1 = 233; % skin panel width, mm
tau_cr = Ks * E * 1e9 * (t2 / b1) ^ 2 * 1e-6; % critical shear buckling stress, MPa
sigma_0 = 430; % Compressive stress of skin, MPa (FROM SKIN)
sigma_cr = 564; % Critical compressive stress, MPa (FROM SKIN)
R_c = sigma_0 / sigma_cr; % Compressive stress ratio
R_s = tau_0 / tau_cr; % Shear stress ratio
val = R_s^2 + R_c; % combined stress ratio

% ribs
M = 43627170; % bending moment, Nm
c = 13.99; % chord, m
t_e = 17.84; % equivalent panel thickness, mm (FROM SKIN)
s = 2; % rib spacing, m
I = (c * (t_e / 1000) ^ 3 / 12 + c * (t_e / 1000) * (b2 / 2) ^ 2); % 2nd moment of area of panel, m^4
F = M^2 * s * b2 * (t_e / 1000) * c / (2 * E * 1e9) / (I ^ 2); % Crush load, N
t_r = 7.3; % design rib thickness, mm
sigma_c = (F / t_r) * (c / 1000); % crush stress, MPa
sigma_b = F / (t_r * c * 1e3); % buckling stress, MPa
t_r2 = F / (sigma_y * c * 1e3); % design rib thickness (again???), MPa