clear
clc

% Aircraft Parameters
l_plane = 77.8; % length of aircraft
W0 = 354000*9.81; % aircraft weight
W_wing = 27100*9.81; % wing weight
W_engine = 6200*9.81; % engine weight
W_fuel = 169*804*9.81; % fuel weight
W_fuselage = 31300*9.81; % fuselage weight
W_tail = (1970+1900)*9.81; % tailplane weight
x_wing = 34.8; % wing location
x_tail = 72.1; % tail location
T_max = 321.6*1000*4; % maximum takeoff thrust of 4 engines
MAC = 7.41; % mean aerodynamic chord
S = 482; % S_ref
x_cg = 34.78; % aircraft cg
C_M = -0.131; % pitching moment

% LOAD CASE 1: 
% Symmetric flight at the Ultimate Load factor, 
% evaluated at manoeuvre VA and dive VD speeds
% VA = 254knots (Cl = 2.8123)
% VD = 333knots (Cl = 1.6341)
V_A = 254*0.514444; % manoeuvre speed
V_D = 333*0.514444; % dive speed
n = 2.5*1.5; % ultimate load factor

%%% 1. Reaction Force from Tail

xw_le = 32; % wing leading edge location
xw_AC = xw_le+0.25*MAC; % aerodynamic centre of wing
l1 = abs(x_cg-xw_AC); 
l2 = 72.1-xw_AC; % distance of tail AC (assume AC coincide with CG) from aircraft AC
[~,a,~,rho,~,~] = atmosisa(10000); 
M_A = 0.5*rho*(V_A)^2*S*MAC*C_M; % pitching moment at VA
M_D = 0.5*rho*(V_D)^2*S*MAC*C_M; % pitching moment at VD
L_tail_A = (n*W0*l1-M_A)/l2; % lift produce by tail at VA (downwards positive)
L_tail_D = (n*W0*l1-M_D)/l2; % lift produce by tail at VD (downwards positive)

%%% 2. Reaction Force from Wing

%%%%%% i. Lift
L_wing = W0*n; 
x_fspar = 33.4; % location of wing front spar
x_rspar = 40.66; % location of wing rear spar
l_F = abs(x_cg-x_fspar);
l_R = abs(x_rspar-x_cg);
% solve two simultaneous equation of moment about cg to find reaction at
% front and rear spar
% VA
A = [l_F, -l_R;
    1, 1];
B = [-L_tail_A+l2;
    n*W0+L_tail_A];
x = A^(-1)*B;
R_F_A = x(1) ; % reaction force on front spar
R_R_A = x(2); % reaction force on rear spar
% VD
A = [l_F, -l_R;
    1, 1];
B = [-L_tail_D+l2;
    n*W0+L_tail_D];
x = A^(-1)*B;
R_F_D = x(1) ; % reaction force on front spar
R_R_D = x(2); % reaction force on rear spar

%%%%%% ii. Weight
W_wing_total = (W_wing+W_engine+W_fuel)*n;

%%%%%% iii. Thrust (neglect)


%%% 3. Fuselage Weight

q_fuselage = W_fuselage/l_plane; % load per unit length 
q_pass = 45*47.88; % load per unit floor area
q_pass = q_pass*6.34; % load per unit fuselage length
q_lugg = 20*157.087; % load per unit cargo volume;


% discretize
n_discrete = 1000;
dl = l_plane/n_discrete;
distance = linspace(0,l_plane,n_discrete);
load1_A = zeros(1000);
load1_D = zeros(1000);

% populate load values (VA)
load1_A(:) = -(q_fuselage+q_pass)*n*dl;
load1_A(round(x_fspar/dl)) = load1_A(round(x_fspar/dl))+ R_F_A-W_wing_total/2;
load1_A(round(x_rspar/dl)) = load1_A(round(x_rspar/dl))+ R_R_A-W_wing_total/2;
load1_A(round(x_tail/dl)) = load1_A(round(x_tail/dl))- L_tail_A-W_tail;

% populate load values (VA)
load1_D(:) = -(q_fuselage+q_pass)*n*dl;
load1_D(round(x_fspar/dl)) = load1_D(round(x_fspar/dl))+ R_F_D-W_wing_total/2;
load1_D(round(x_rspar/dl)) = load1_D(round(x_rspar/dl))+ R_R_D-W_wing_total/2;
load1_D(round(x_tail/dl)) = load1_D(round(x_tail/dl))- L_tail_D-W_tail;

% Plot Load Distribution
figure(1)
plot(distance,load1_A,'LineStyle',':');
hold on
plot(distance,load1_D,'LineStyle','-');
title("Weight distribution along fuselage");
xlabel("Length (m)");
ylabel("Weight distribution (N)");
grid on


% LOAD CASE 2: OEI
% 1.885e5 N act at CG of Vertical Stabliser
% create torsion on the fuselage

% LOAD CASE 3: Landing with main gears only
% assume no lift on wing and tailplane
% Landing Load Factor
Vs = 12; % sink speed 12ft/sec
% Tire spring constant = 160000 lb/ft
% Tire deflection 20.7 inch
% L = 2/3 W
% struct stroke 0.25m
% strut efficiency, 0.8 ot 0.85


%%
% calculate shear
shear = zeros(n_discrete);
for i = 2:n_discrete
    shear(i) = shear(i-1) - load1_A(i);
end

% Plot Shear Distribution
figure(2)
plot(distance,shear);
title("Shear distribution along fuselage");
xlabel("Length (m)");
ylabel("Shear distribution (N)");
grid on

% calculate moment
moment = zeros(n_discrete);
for i = 2:n_discrete
    moment(i) = moment(i-1) + shear(i);
end

% Plot Shear Distribution
figure(3)
plot(distance,moment);
title("Moment distribution along fuselage");
xlabel("Length (m)");
ylabel("Moment distribution (N)");
grid on






