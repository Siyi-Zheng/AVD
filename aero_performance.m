clear
clc

AR = 8.77;          %aspect ratio
Cla = 6.74;         %lift curve slope of aerofoil
d = 6.34;           %diameter of fuselage
b = 65;             %wing span
S_exp = 440;        %exposed area
S_ref = 482;        %reference area
sweep = 26;         % sweep at max thickness pt.
le_sweep = 30;      %leading edge sweep

CLa1 = [];
CLa2 = [];
CLa3 = [];
M_list1 = [];
M_list2 = [];
for M = 0:0.01:2

    beta = (1 - M ^ 2) ^ 0.5;
    eta = Cla / (2 * pi); % remove the beta term as Cla has it already
    F = 1.07 * (1 + d / b) ^ 2;

    CLa_subsonic3 = 2 * pi * AR * S_exp * F / (2 + (4 + (AR * beta / eta)...
        ^ 2 *(1 + (tand(sweep) / beta) ^ 2)) ^ 0.5) / S_ref;
    CLa_subsonic = 2*pi/sqrt(1-M^2);
    CLa_supersonic = 4/sqrt(-1+M^2);

    if M < 1
        CLa1 = [CLa1 CLa_subsonic];
        CLa3 = [CLa3 CLa_subsonic3];
        M_list1 = [M_list1 M];
    elseif M > 1 / cosd(sweep)
        CLa2 = [CLa2, CLa_supersonic];
        M_list2 = [M_list2 M];
    end
end

figure
plot(M_list1, CLa1,"red", LineWidth=1.5);
hold on
plot(M_list2, CLa2,"blue", LineWidth=1.5);
plot(M_list1, CLa3,"magenta", LineWidth=1.5);
xline(0.83, LineWidth=1.5);
xlabel("M");
ylabel("CLa")
ylim([0 12]);
legend("Subsonic Theoretical","Supersonic Theoretical","Subsonic Raymer","Subsonic 2");
grid on




% drag prediction IN CRUISE

%=========================================================================

% CD0
% for smooth paint, k = 6.34e-6 m - so we are below the cutoff

% fuselage
Re = 434e6;
M = 0.83;
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

CL_cruise = 0.53;
AR_h = 5.88;
e = 0.85; %oswald efficiency factor
K = 1 / (pi * AR * e);
K_h = 1 / (pi * AR_h * e);
eta_h = 0.9; % tailplane efficiency
CL_h = 0; % we need stability analysis to get this one
CDi = K * CL_cruise ^ 2 + eta_h * K_h * CL_h ^ 2;

% also get the change in CDi from flaps

CDi_flaps = 0; % idk what k_h is so we can sort this out later




%=========================================================================

%IMPORTANT VALUES!!!!!
%CLa = 7.738 FOR THE WING FOR M = 0.83
