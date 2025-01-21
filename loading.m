clear all
clc

% --- CRUISE ---

% symmetric loading
rho = 1.227; % sea level density (kg/m^3)
mac = 7.4; % mean aero chord (m)
Cla = 6.29; % lift curve slope (1/rad)
wing_loading = 353385 / 482 * 9.81; % wing loading (kg/m^2)
n_max = 2.5; % max. load factor
n_min = -1.0; % min. load factor
v_max = 200; % max. speed (m/s) --- how do i get this properly???
CL_max = 1.72; % with flaps
CL_min = -1.5; % idk what this actually is
v_list = 0:1:v_max*1.4; % for plotting the envelope
v_stall = sqrt(2 * wing_loading / (rho * CL_max));

% aero loads
n_max_aero = rho * CL_max / (2 * wing_loading) * v_list .^ 2;
n_min_aero = rho * CL_min / (2 * wing_loading) * v_list .^ 2;

% gust loads
mu = 2 * wing_loading / (rho * 9.81 * mac * Cla);
U_de_c = 15.2; % gust speed (m/s) at cruise
U_c = U_de_c * 0.88 * mu / (5.3 + mu); % actual gust velocity
n_max_gust_c = 1 + rho * U_c * v_list * Cla / (2 * wing_loading);
n_min_gust_c = 1 - rho * U_c * v_list * Cla / (2 * wing_loading);
U_de_d = 7.6; % gust speed (m/s) in a dive
U_d = U_de_d * 0.88 * mu / (5.3 + mu); % actual gust velocity
n_max_gust_d = 1 + rho * U_d * v_list * Cla / (2 * wing_loading);
n_min_gust_d = 1 - rho * U_d * v_list * Cla / (2 * wing_loading);

% plotting
% aero loads
plot(v_list, n_max_aero, color="black", linestyle="--")
hold on
plot(v_list, n_min_aero, color="black", linestyle="--")
yline(n_max, color="black")
yline(n_min, color="black")
xline(v_max, color="black")
xline(v_stall, color="black")
% gust loads
plot(v_list, n_max_gust_c, color="black", linestyle=":")
plot(v_list, n_min_gust_c, color="black", linestyle=":")
plot(v_list, n_max_gust_d, color="black", linestyle=":")
plot(v_list, n_min_gust_d, color="black", linestyle=":")
% formatting
xlim([min(v_list) max(v_list)])
ylim([-1.5 3])
xlabel("V_{EAS} (m/s)")
ylabel("Load factor")
grid on