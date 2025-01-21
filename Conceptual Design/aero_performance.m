clear
clc


% wave drag
M_crit = 0.88;
LF_DD = 0.875;
CL_des = 0.54;

M_dd_0 = M_crit + 0.08;
M_dd = M_dd_0 * LF_DD - 0.05 * CL_des;


AR = 8.77;          %aspect ratio
Cla = 6.74;         %lift curve slope of aerofoil
d = 6.34;           %diameter of fuselage
b = 65;             %wing span
S_exp = 440;        %exposed area
S_ref = 482;        %reference area
sweep = 26;         % sweep at max thickness pt.
le_sweep = 30;      %leading edge sweep (degrees)

CLa2 = [];
CLa3 = [];
M_list1 = [];
M_list2 = [];

for M = 0:0.01:1.3

    beta = (1 - M ^ 2) ^ 0.5;
    eta = Cla / (2 * pi); % remove the beta term as Cla has it already
    % if eta > 0.98
    %     eta = 0.98; % correction from Raymer
    % end
    F = 1.07 * (1 + d / b) ^ 2;
    if F * S_exp / S_ref > 1
        S_ref = F * S_exp / 0.98;
    end

    CLa_subsonic3 = 2 * pi * AR * S_exp * F * cosd(le_sweep) / (2 + (4 + (AR * beta / eta)...
        ^ 2 *(1 + (tand(sweep) / beta) ^ 2)) ^ 0.5) / S_ref;
    CLa_supersonic = 4/sqrt(-1+M^2);
    CLa_sonic = 2 * pi * AR * S_exp * F * cosd(le_sweep) / (2 + (4 + (AR * 1e-6 / eta)...
        ^ 2 *(1 + (tand(sweep) / 1e-6) ^ 2)) ^ 0.5) / S_ref;

    if M < 1 / cosd(sweep) && M ~= 1
        CLa3 = [CLa3 CLa_subsonic3];
        M_list1 = [M_list1 M];
    elseif M == 1
        CLa3 = [CLa3 CLa_sonic];
        M_list1 = [M_list1 M];
    end
    if M > 1 / cosd(sweep) - 0.1
        CLa2 = [CLa2, CLa_supersonic];
        M_list2 = [M_list2 M];
    end
end

% plot the results
figure
hold on
plot(M_list2, CLa2,"blue", LineWidth=1.5);
plot(M_list1, CLa3,"red", LineWidth=1.5);
xregion(0.85, 1.2)
xline(0.83, label="Cruise");
xline(0.23, label="Takeoff");
xline(0.27, label="Landing");
ylim([3 10])
xlabel("M");
ylabel("C_{L_α}")
xlim([0 1.3]);
legend("Supersonic","Subsonic","Transonic Regime");
grid on




% drag prediction IN CRUISE

%=========================================================================

% CD0
% for smooth paint, k = 6.34e-6 m - so we are below the cutoff

% fuselage
Re = 434e6;
M = 0.83;
S_ref = 482;
f = 80/6.34; % fineness ratio
Cf_fus = 0.455 / (log10(Re) ^ 2.58 * (1 + 0.144 * M^2) ^ 0.65);
FF_fus = 0.9 + 5 * f ^ -1.5 + f / 400; % form factor
S_wet_fus = 1.38e3;
CD0_fus = S_wet_fus / S_ref * Cf_fus * FF_fus;

% wings
Re = 40.1e6;
M = 0.83;
xc_max = 0.36; % x/c of max. thickness
tc_max = 0.14; % max thickness-to-chord ratio
Cf_w = 0.455 / (log10(Re) ^ 2.58 * (1 + 0.144 * M^2) ^ 0.65);
FF_w = (1 + 0.6/xc_max * tc_max + 100 * (tc_max) ^ 4) * ...
    (1.34 * M ^ 0.18 * cosd(sweep) ^ 0.28); % form factor
S_wet_w = S_exp * 2;
CD0_w = S_wet_w / S_ref * Cf_w * FF_w;

% hstab
Re = 22.8e6;
M = 0.83;
S_h = 51.7;
xc_max = 0.298; % x/c of max. thickness
tc_max = 0.12; % max thickness-to-chord ratio
Cf_h = 0.455 / (log10(Re) ^ 2.58 * (1 + 0.144 * M^2) ^ 0.65);
FF_h = (1 + 0.6/xc_max * tc_max + 100 * (tc_max) ^ 4) * ...
    (1.34 * M ^ 0.18 * cosd(sweep) ^ 0.28); % form factor
S_wet_h = S_h * 2.1;
CD0_h = S_wet_h / S_ref * Cf_h * FF_h;

% vstab
Re = 32.9e6;
M = 0.83;
S_v = 27.7;
xc_max = 0.3; % x/c of max. thickness
tc_max = 0.12; % max thickness-to-chord ratio
Cf_v = 0.455 / (log10(Re) ^ 2.58 * (1 + 0.144 * M^2) ^ 0.65);
FF_v = (1 + 0.6/xc_max * tc_max + 100 * (tc_max) ^ 4) * ...
    (1.34 * M ^ 0.18 * cosd(sweep) ^ 0.28); % form factor
S_wet_v = S_v * 2.1;
CD0_v = S_wet_v / S_ref * Cf_v * FF_v;

% nacelles
Re = 35.2e6;
M = 0.83;
f = 4.64; % fineness ratio
Q = 1.3; % nacelle-wing interference
S_wet_nac = 69.5 * 4; % approximate as a cylinder, and there are 4
Cf_nac = 0.455 / (log10(Re) ^ 2.58 * (1 + 0.144 * M^2) ^ 0.65);
FF_nac = 1 + 0.35 / f; % form factor
CD0_nac = S_wet_nac / S_ref * Cf_nac * FF_nac * Q;

% fuselage upsweep
A_fus = 31.6; % cross-sectional area
beta = 7 * pi/180; % upsweep angle (rad)
CD_upsweep = 3.83 * A_fus / S_ref * beta ^ 2.5;

% get clean CD0 (no gear or flaps)
CD0_clean = CD0_fus + CD0_w + CD0_h + CD0_v + CD0_nac + CD_upsweep;

% get CD changes from different configs:
% landing gear
A_uc_front = 13.9; % landing gear area from front - get a better value?
CD_uc = 2.25 * A_uc_front / S_ref;

% flaps
flap_extent = 0.65;
CD_flaps_TO = 0.0023 * flap_extent * 15; % 15deg for takeoff
CD_flaps_land = 0.0024 * flap_extent * 40; % 40deg for landing

% CDi
% we need to finish stability to get all these values i think

CL_cruise = 2.85;
AR_h = 5.88;
e = 0.85; %oswald efficiency factor
K = 1 / (pi * AR * e);
K_h = 1 / (pi * AR_h * e);
eta_h = 0.9; % tailplane efficiency
CL_h = -0.27; % we need stability analysis to get this one
CDi = K * CL_cruise ^ 2 + eta_h * K_h * CL_h ^ 2;

% also get the change in CDi from flaps

kf = 1.4;
delta_CL = 0.139 * 2.5;

CDi_flaps = kf ^ 2 * delta_CL ^ 2 + cosd(sweep); % idk what k_f is so we can sort this out later

%=========================================================================

% drag prediction AT TAKEOFF

%=========================================================================

% CD0
% for smooth paint, k = 6.34e-6 m - so we are below the cutoff

% fuselage
Re = 423e6;
M = 0.23;
f = 80/6.34; % fineness ratio
Cf_fus = 0.455 / (log10(Re) ^ 2.58 * (1 + 0.144 * M^2) ^ 0.65);
FF_fus = 0.9 + 5 * f ^ -1.5 + f / 400; % form factor
S_wet_fus = 1.38e3;
CD0_fus_TO = S_wet_fus / S_ref * Cf_fus * FF_fus;

% wings
Re = 39.1e6;
M = 0.23;
xc_max = 0.36; % x/c of max. thickness
tc_max = 0.14; % max thickness-to-chord ratio
Cf_w = 0.455 / (log10(Re) ^ 2.58 * (1 + 0.144 * M^2) ^ 0.65);
FF_w = (1 + 0.6/xc_max * tc_max + 100 * (tc_max) ^ 4) * ...
    (1.34 * M ^ 0.18 * cosd(sweep) ^ 0.28); % form factor
S_wet_w = S_exp * 2;
CD0_w_TO = S_wet_w / S_ref * Cf_w * FF_w;

% hstab
Re = 22.2e6;
M = 0.23;
S_h = 51.7;
xc_max = 0.298; % x/c of max. thickness
tc_max = 0.12; % max thickness-to-chord ratio
Cf_h = 0.455 / (log10(Re) ^ 2.58 * (1 + 0.144 * M^2) ^ 0.65);
FF_h = (1 + 0.6/xc_max * tc_max + 100 * (tc_max) ^ 4) * ...
    (1.34 * M ^ 0.18 * cosd(sweep) ^ 0.28); % form factor
S_wet_h = S_h * 2.1;
CD0_h_TO = S_wet_h / S_ref * Cf_h * FF_h;

% vstab
Re = 32e6;
M = 0.23;
S_v = 27.7;
xc_max = 0.3; % x/c of max. thickness
tc_max = 0.12; % max thickness-to-chord ratio
Cf_v = 0.455 / (log10(Re) ^ 2.58 * (1 + 0.144 * M^2) ^ 0.65);
FF_v = (1 + 0.6/xc_max * tc_max + 100 * (tc_max) ^ 4) * ...
    (1.34 * M ^ 0.18 * cosd(sweep) ^ 0.28); % form factor
S_wet_v = S_v * 2.1;
CD0_v_TO = S_wet_v / S_ref * Cf_v * FF_v;

% nacelles
Re = 34.3e6;
M = 0.23;
f = 4.64; % fineness ratio
Q = 1.3; % nacelle-wing interference
S_wet_nac = 69.5 * 4; % approximate as a cylinder, and there are 4
Cf_nac = 0.455 / (log10(Re) ^ 2.58 * (1 + 0.144 * M^2) ^ 0.65);
FF_nac = 1 + 0.35 / f; % form factor
CD0_nac_TO = S_wet_nac / S_ref * Cf_nac * FF_nac * Q;

% get takeoff CD0
CD0_TO = CD0_fus_TO + CD0_w_TO + CD0_h_TO + CD0_v_TO + CD0_nac_TO +...
    CD_upsweep + CD_flaps_TO + CD_uc;

%=========================================================================

% drag prediction AT LANDING

%=========================================================================

% CD0
% for smooth paint, k = 6.34e-6 m - so we are below the cutoff

% fuselage
Re = 500e6;
M = 0.27;
f = 80/6.34; % fineness ratio
Cf_fus = 0.455 / (log10(Re) ^ 2.58 * (1 + 0.144 * M^2) ^ 0.65);
FF_fus = 0.9 + 5 * f ^ -1.5 + f / 400; % form factor
S_wet_fus = 1.38e3;
CD0_fus_land = S_wet_fus / S_ref * Cf_fus * FF_fus;

% wings
Re = 46.2e6;
M = 0.27;
xc_max = 0.36; % x/c of max. thickness
tc_max = 0.14; % max thickness-to-chord ratio
Cf_w = 0.455 / (log10(Re) ^ 2.58 * (1 + 0.144 * M^2) ^ 0.65);
FF_w = (1 + 0.6/xc_max * tc_max + 100 * (tc_max) ^ 4) * ...
    (1.34 * M ^ 0.18 * cosd(sweep) ^ 0.28); % form factor
S_wet_w = S_exp * 2;
CD0_w_land = S_wet_w / S_ref * Cf_w * FF_w;

% hstab
Re = 26.2e6;
M = 0.27;
S_h = 51.7;
xc_max = 0.298; % x/c of max. thickness
tc_max = 0.12; % max thickness-to-chord ratio
Cf_h = 0.455 / (log10(Re) ^ 2.58 * (1 + 0.144 * M^2) ^ 0.65);
FF_h = (1 + 0.6/xc_max * tc_max + 100 * (tc_max) ^ 4) * ...
    (1.34 * M ^ 0.18 * cosd(sweep) ^ 0.28); % form factor
S_wet_h = S_h * 2.1;
CD0_h_land = S_wet_h / S_ref * Cf_h * FF_h;

% vstab
Re = 37.8e6;
M = 0.27;
S_v = 27.7;
xc_max = 0.3; % x/c of max. thickness
tc_max = 0.12; % max thickness-to-chord ratio
Cf_v = 0.455 / (log10(Re) ^ 2.58 * (1 + 0.144 * M^2) ^ 0.65);
FF_v = (1 + 0.6/xc_max * tc_max + 100 * (tc_max) ^ 4) * ...
    (1.34 * M ^ 0.18 * cosd(sweep) ^ 0.28); % form factor
S_wet_v = S_v * 2.1;
CD0_v_land = S_wet_v / S_ref * Cf_v * FF_v;

% nacelles
Re = 40.5e6;
M = 0.27;
f = 4.64; % fineness ratio
Q = 1.3; % nacelle-wing interference
S_wet_nac = 69.5 * 4; % approximate as a cylinder, and there are 4
Cf_nac = 0.455 / (log10(Re) ^ 2.58 * (1 + 0.144 * M^2) ^ 0.65);
FF_nac = 1 + 0.35 / f; % form factor
CD0_nac_land = S_wet_nac / S_ref * Cf_nac * FF_nac * Q;

% get takeoff CD0
CD0_land = CD0_fus_land + CD0_w_land + CD0_h_land + CD0_v_land +...
    CD0_nac_land + CD_upsweep + CD_flaps_land + CD_uc;



%IMPORTANT VALUES!!!!!
% CLa = 5.581 for cruise
% CLa = 4.233 for landing
% CLa = 4.205 for takeoff
% CD0 clean = 0.0161
% CD0 takeoff = 0.1026
% CD0 landing = 0.1424
