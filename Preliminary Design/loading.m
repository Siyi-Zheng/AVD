clear all
clc

% symmetric loading
rho = 1.227; % sea level density (kg/m^3)
mac = 7.4; % mean aero chord (m)
Cla = 6.29; % lift curve slope (1/rad)
wing_loading = 353385 / 482 * 9.81; % wing loading (kg/m^2)
n_max = 2.5; % max. load factor
n_min = -1.0; % min. load factor
v_c = 137; % cruise speed (m/s EAS)
v_d = v_c / 0.8; % dive speed (m/s) using FAR 25
CL_max = 1.72;
CL_max_flaps = 2.85;
CL_min = -1.01; % idk what this is actually
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
U_ref_c = 11.06;
F_gz = 1 - 43000/250000;
R1 = 0.85;
R2 = (353385 - 176090) / 353385;
F_gm = sqrt(R2 * tan(pi * R1 / 4));
F_g = 0.5 * (F_gz + F_gm);
U_ref_d = U_ref_c * 0.5;
H = 107;
U_c = U_ref_c * F_g * (H / 107) ^ (1/6);
n_max_gust_c = 1 + (rho * U_c * v_list * Cla / (2 * wing_loading));
n_min_gust_c = 1 - (rho * U_c * v_list * Cla / (2 * wing_loading));
U_d = U_ref_d * F_g * (H / 107) ^ (1/6);
n_max_gust_d = 1 + (rho * U_d * v_list * Cla / (2 * wing_loading));
n_min_gust_d = 1 - (rho * U_d * v_list * Cla / (2 * wing_loading));

% plotting
% aero loads
plot(v_list*1.944, n_max_aero, color="black", linestyle="--")
hold on
plot(v_list*1.944, n_min_aero, color="black", linestyle="--")
yline(n_max, color="black", linestyle="--")
yline(n_min, color="black", linestyle="--")
xline(v_a*1.944, color="black", linestyle="--")
xline(v_c*1.944, color="black", linestyle="--")
xline(v_d*1.944, color="black", linestyle="--")
xline(v_stall*1.944, color="black", linestyle="--")
xline(v_stall_flaps*1.944, color="black", linestyle="--")
xline(v_b_d*1.944, color="black", linestyle="--")
xline(v_flaps*1.944, color="black", linestyle="--")

% gust loads
plot(v_list*1.944, n_max_gust_c, color="black", linestyle=":")
plot(v_list*1.944, n_min_gust_c, color="black", linestyle=":")
plot(v_list*1.944, n_max_gust_d, color="black", linestyle=":")
plot(v_list*1.944, n_min_gust_d, color="black", linestyle=":")

% final envelope
plot(v_list_flaps*1.944, n_max_flaps, color="black", linewidth=2)
plot([v_stop_flaps*1.944 v_flaps*1.944], [2 2], color="black", linewidth=2)
n_flaps_max = rho * CL_max / (2 * wing_loading) * v_flaps .^ 2;
plot([v_flaps*1.944 v_flaps*1.944], [n_flaps_max 2], color="black", linewidth=2)
v_aero_str = sqrt(5 * wing_loading / (rho * CL_max));
v_aero_str2 = sqrt(-2 * wing_loading / (rho * CL_min));
v_aero_list = linspace(v_flaps, v_aero_str, 100);
n_aero_list = rho * CL_max / (2 * wing_loading) * v_aero_list .^ 2;
plot(v_aero_list*1.944, n_aero_list, color="black", linewidth=2)
plot([v_aero_str*1.944 v_d*1.944], [2.5 2.5], color="black", linewidth=2)
plot([v_d*1.944 v_d*1.944], [2.5 0], color="black", linewidth=2)
plot([v_aero_str2*1.944 v_c*1.944], [-1 -1], color="black", linewidth=2)
v_aero_list2 = linspace(0, v_aero_str2, 100);
n_aero_list2 = rho * CL_min / (2 * wing_loading) * v_aero_list2 .^ 2;
plot(v_aero_list2*1.944, n_aero_list2, color="black", linewidth=2)
v_diagonal_list = linspace(v_c, v_d, 100);
n_diagonal_list = linspace(-1, 0, 100);
plot(v_diagonal_list*1.944, n_diagonal_list, color="black", linewidth=2)

% formatting
xlim([min(v_list)*1.944 max(v_list)*1.944])
ylim([-1.5 3.5])
xlabel("Airspeed (KEAS)")
ylabel("Load factor")
grid on