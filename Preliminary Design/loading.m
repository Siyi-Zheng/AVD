clear all
clc

% symmetric loading
rho = 1.227; % sea level density (kg/m^3)
mac = 7.4; % mean aero chord (m)
Cla = 6.29; % lift curve slope (1/rad)
wing_loading = 353385 / 482 * 9.81; % wing loading (kg/m^2)
n_max = 2.5; % max. load factor
n_min = -1.0; % min. load factor
v_c = 245; % cruise speed (m/s EAS)
v_d = v_c * (0.9 / 0.83); % dive speed (m/s) using FAR 25
CL_max = 1.72;
CL_max_flaps = 2.85;
CL_min = -1; % idk what this is actually
v_list = linspace(0, v_d+20, 100); % for plotting the envelope
v_stall = sqrt(2 * wing_loading / (rho * CL_max));
v_stall_flaps = sqrt(2 * wing_loading / (rho * CL_max_flaps));
v_stop_flaps = sqrt(4 * wing_loading / (rho * CL_max_flaps));
v_flaps = 1.6 * v_stall_flaps;
v_list_flaps = linspace(0, v_stop_flaps, 100);
v_a = v_stall * sqrt(2.5);

% aero loads
n_max_aero = rho * CL_max / (2 * wing_loading) * v_list .^ 2;
n_min_aero = rho * CL_min / (2 * wing_loading) * v_list .^ 2;
n_max_flaps = rho * CL_max_flaps / (2 * wing_loading) * v_list_flaps .^ 2;

% V_b
v_b_d = v_stall * (1 + (rho * 20 * v_stall * Cla) / (2 * wing_loading)) ^ 0.5;

% gust loads
mu = 2 * wing_loading / (rho * 9.81 * mac * Cla);
U_de_c = 15.2; % gust speed (m/s) at cruise
U_c = U_de_c * 0.88 * mu / (5.3 + mu); % actual gust velocity
n_max_gust_c = 1 + (rho * U_c * v_list * Cla / (2 * wing_loading));
n_min_gust_c = 1 - (rho * U_c * v_list * Cla / (2 * wing_loading));
U_de_d = 7.6; % gust speed (m/s) in a dive
U_d = U_de_d * 0.88 * mu / (5.3 + mu); % actual gust velocity
n_max_gust_d = 1 + (rho * U_d * v_list * Cla / (2 * wing_loading));
n_min_gust_d = 1 - (rho * U_d * v_list * Cla / (2 * wing_loading));

% plotting
% aero loads
plot(v_list, n_max_aero, color="black", linestyle="--")
hold on
plot(v_list, n_min_aero, color="black", linestyle="--")
yline(n_max, color="black", linestyle="--")
yline(n_min, color="black", linestyle="--")
xline(v_a, color="black", linestyle="--", label="V_A")
xline(v_c, color="black", linestyle="--", label="V_C")
xline(v_d, color="black", linestyle="--", label="V_D")
xline(v_stall, color="black", linestyle="--", label="V_{S1}")

% gust loads
plot(v_list, n_max_gust_c, color="black", linestyle=":")
plot(v_list, n_min_gust_c, color="black", linestyle=":")
plot(v_list, n_max_gust_d, color="black", linestyle=":")
plot(v_list, n_min_gust_d, color="black", linestyle=":")

% final envelope
plot(v_list_flaps, n_max_flaps, color="black", linewidth=2)
plot([v_stop_flaps v_flaps], [2 2], color="black", linewidth=2)
n_flaps_max = rho * CL_max / (2 * wing_loading) * v_flaps .^ 2;
plot([v_flaps v_flaps], [n_flaps_max 2], color="black", linewidth=2)
v_aero_str = sqrt(5 * wing_loading / (rho * CL_max));
v_aero_str2 = sqrt(-2 * wing_loading / (rho * CL_min));
v_aero_list = linspace(v_flaps, v_aero_str, 100);
n_aero_list = rho * CL_max / (2 * wing_loading) * v_aero_list .^ 2;
plot(v_aero_list, n_aero_list, color="black", linewidth=2)
plot([v_aero_str v_d], [2.5 2.5], color="black", linewidth=2)
plot([v_d v_d], [2.5 -1], color="black", linewidth=2)
plot([v_aero_str2 v_d], [-1 -1], color="black", linewidth=2)
v_aero_list2 = linspace(0, v_aero_str2, 100);
n_aero_list2 = rho * CL_min / (2 * wing_loading) * v_aero_list2 .^ 2;
plot(v_aero_list2, n_aero_list2, color="black", linewidth=2)

% formatting
xlim([min(v_list) max(v_list)])
ylim([-1.5 3.5])
xlabel("V_{EAS} (m/s)")
ylabel("Load factor")
grid on