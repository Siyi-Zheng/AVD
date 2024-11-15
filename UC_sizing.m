clc
clear

%Length
L_plane = 80.6; % length of plane in METRE
L_ws = 65; % wing span in METRE

%Weights
W_o = 390057;        %MTOW in KG
W_ng = 0.085*W_o;      %Nose gear weight carrying capability
W_mg = W_o - W_ng;   %Main gear weight carrying capability

MAC = 7.41; %mean aerodynamic chord in METRE

Wo_lbs = W_o*2.20462; %Wo in lbs (errikos just being extremely annoying with values)
%This value is 85 9927 lbs

%Moments
x_wing = 37; % distance of wing leading edge in METRE
x_cg = 37 + 0.4*MAC;                        % distance of cg (estimate)
x_mg = 43;                                %x distance of main gears (estimate)
x_ng = (W_o*x_cg - W_mg*x_mg)/W_ng;         %x distance of nose gears


%Main Gear Height
beta = 15; % tipback angle in degree
gamma = 15; % rotation clearance angle
H = tan(deg2rad(gamma))*(L_plane-x_mg-200/39.37); % correction for the rear main gear
tip_back_angle = rad2deg(tan((x_mg-x_cg)/H)); % between 15 ~ 25 degree
%check wing tip do not hit the ground at 5 degree roll and AOA 90% Cl_max


%Main Gear Separation
L_mg = 220 * 0.0254 ; % horizontal distance from main landing gear to plane centre axis in METRES
% angle between line connecting main and nose gear to plane centre axis
a = rad2deg(atan(L_mg/(x_mg-x_ng)));
d = (x_cg-x_ng)*sin(deg2rad(a)); % static ground line
a_overturn = rad2deg(atan(H/d)); % should be below 63 degree


% Loading (IN KG)
% Assume aftmost CG is x_cg
% Assume foremost CG is 20~25% MAC in front of aftmost CG
x_cg_foremost = x_cg - 0.20*MAC; % foremost cg location in METRE
W_mg_max = W_o*(x_cg-x_ng)/(x_mg-x_ng) * 1.07; % max static loading on main gear in KG
W_ng_max = W_o*(x_mg-x_cg_foremost)/(x_mg-x_ng) * 1.07; % max static loading on nose gear in KG
H_feet = H*3.28084; % METRE to FEET
W_ng_braking_lbs = 10*H_feet*Wo_lbs/(32.2*(x_mg-x_ng) * 3.28084) * 1.07; % dynamic braking load in lbs
W_ng_braking = 0.453592*W_ng_braking_lbs; % lbs to KG

%Nose Gear: one struct, two wheels
%Tire Sizing (statistical approach): Nose wheel
W_ng_lbs = W_ng_max*2.20462; % convert loading from KG to lbs
W_nw_lbs = (W_ng_lbs+W_ng_braking_lbs)/2; %Weight on each nose wheel, factor 1.07 from errikos's slide
D_nw = 1.63*W_nw_lbs^0.315; %Diameter of nose wheel in INCH
Wth_nw = 0.1043*W_nw_lbs^0.48; %Width of nose wheel in INCH

%Main Gear: four struct, each with 4 wheels (twin tandem)
%Tire Sizing (statistical approach): main wheel
W_mg_lbs = W_mg*2.20462; % convert loading from KG to lbs
W_mw_lbs = W_mg_lbs*1.07/16; %Weight on each main wheel in lbs
D_mw = 1.63*W_mw_lbs^0.315; %Diameter of main wheel in INCH
Wth_mw =  0.1043*W_mw_lbs^0.48; %Width of main wheel in INCH


%Tire Selection (From a table on Raymer's book), based on Max loading
% Main Wheel
%Three-part Name Tire: 52 x 20.5-23
R_r_m = 21.3; % Rolling Radius in INCH
w_tire_m = 20.5;
d_tire_m = 52;

% Nose Wheel
R_r_n = 21.3; % Rolling Radius in INCH
w_tire_n = 20.5; % width
d_tire_n = 52; % diameter


%Tire Pressure
% Main Tire
A_p_mw = 2.3*sqrt(w_tire_m*d_tire_m)*(d_tire_m/2 - R_r_m); %wheel contact area in inch^2
P_mw = W_mw_lbs/A_p_mw; % tire pressure in PSI

% Nose Tire
A_p_nw = 2.3*sqrt(w_tire_n*d_tire_n)*(d_tire_n/2 - R_r_n); %wheel contact area in inch^2
P_nw = W_nw_lbs/A_p_nw; % tire pressure in PSI


%Brake
W_landing = 0.5*W_o; % account for emergency landing short aft takeoff in KG
KE = 0.5*W_landing * 70.4^2; % total KE required to absorb by the brake
KE_per_wheel = KE/16; % KE required for each wheel with break


%Stroke
V_vertical = 12; % vertical touchdown speed in feet per second
N_gear = 2.7; % gear loading factor 2.7 - 3
S_T = D_mw/2 - R_r_m; % stroke of tire (half radius - rolling radius)
n_tire = 0.47; % shock efficiency of tire
n = 0.85; % shock absorber efficiency
% unit conversion
V_vertical = V_vertical*0.3048; % feet to metre
S_T = S_T*0.0254; % INCH to metre
S = (V_vertical)^2 / (2*9.81*n*N_gear) - n_tire/n * S_T; % should be around 25-30cm


% Oleo Sizing
L_oleo = 2.5*S; % length of oleo
P_oleo = 12415e3; % internal pressure oleo as in Raymers
% Main wheel
W_moleo = W_mg_max*9.81/4; % Load on main wheel in Newton
D_moleo = 1.3 * sqrt((4*W_moleo)/(P_oleo*pi));
D_moleo = 0.04*sqrt(W_moleo*0.224809)*0.0254;
% Nose wheel
W_noleo = (W_ng_max+W_ng_braking)*9.81; % Load on nose wheel in Newton
D_noleo = 0.04*sqrt(W_noleo*0.224809)*0.0254;


%ACN - Aircraft Classification Number
%The follwoing values are generated from the COMFAA software

