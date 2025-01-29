clear
clc

% Aircraft Parameters
l_plane = 77.8; % length of aircraft
W0 = 354000*9.81; % aircraft weight
W_wing = 27100*9.81; % wing weight
W_engine = (26200+5360)*9.81; % engine weight
W_fuel = 76.8*804*9.81; % fuel weight
W_fuselage = 31300*9.81; % fuselage weight
W_htail = 1970*9.81; % horizontal tailplane weight
W_vtail = 1900*9.81; % vertical tailplane weight
W_mlg = 15100/4*9.81; % one main landing gear weight
W_nlg = 550*9.81; % nose landing gear weight
x_wing = 34.8; % wing location
x_htail = 72.1; % horizontal tail location
x_vtail = 71.1; % vertical tail location
x_mlg = 39; % main landing gear location
x_nlg = 7; % nose landing gear location
T_max = 321.6*1000*4; % maximum takeoff thrust of 4 engines
MAC = 7.41; % mean aerodynamic chord
S = 482; % S_ref
x_cg = 35; % aircraft cg
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

xw_AC = 34.5; % aerodynamic centre of wing
l1 = abs(x_cg-xw_AC); 
l2 = 72.1-xw_AC; % distance of tail AC (assume AC coincide with CG) from aircraft AC
[~,a,~,rho,~,~] = atmosisa(0);  
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
B = [-(L_tail_A*l2);
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
% wing empty weight, fuel weight, engine weight, landing gear weight
W_wing_total = (W_wing+W_engine+W_fuel+W_mlg*2);

%%%%%% iii. Thrust (neglect)


%%% 3. Fuselage Weight

q_fuselage = W_fuselage/l_plane; % load per unit length 
q_pass = 45*47.88; % load per unit floor area
q_pass = q_pass*6.34; % load per unit fuselage length
q_lugg = 9428*9.81/l_plane; % load due to cargo ;
q_fuel = 6564.66; % load due to fuselage tank
% front 23.9
% back 32.6
% fuel density 804

% discretize
n_discrete = 1000;
dl = l_plane/n_discrete;
distance = linspace(0,l_plane,n_discrete);
load1_A = zeros(1000,1);
load1_D = zeros(1000,1);

% populate load values (VA)
% uniform load due to fuselage empty weight, luggage, passenger
load1_A(:) = load1_A(:) -(q_fuselage+q_pass+q_lugg)*n*dl; 
% load due to fuel tanl
load1_A(round(23.9/dl):round(32.6/dl)) = load1_A(round(23.9/dl):round(32.6/dl)) - q_fuel*n*dl;
load1_A(round(x_fspar/dl)) = load1_A(round(x_fspar/dl))+ R_F_A-3/4*n*W_wing_total;
load1_A(round(x_rspar/dl)) = load1_A(round(x_rspar/dl))+ R_R_A-1/4*n*W_wing_total;
load1_A(round(x_htail/dl)) = load1_A(round(x_htail/dl)) - L_tail_A - n*W_htail;
load1_A(round(x_vtail/dl)) = load1_A(round(x_vtail/dl)) - n*W_vtail;
load1_A(round(x_mlg/dl)) = load1_A(round(x_mlg/dl)) - n*W_mlg*2;
load1_A(round(x_nlg/dl)) = load1_A(round(x_nlg/dl)) - n*W_nlg;

% populate load values (VD)
% uniform load due to fuselage empty weight, luggage, passenger
load1_D(:) = load1_D(:) -(q_fuselage+q_pass+q_lugg)*n*dl; 
% load due to fuel tanl
load1_D(round(23.9/dl):round(32.6/dl)) = load1_A(round(23.9/dl):round(32.6/dl)) - q_fuel*n*dl;
load1_D(round(x_fspar/dl)) = load1_D(round(x_fspar/dl))+ R_F_D-3/4*n*W_wing_total;
load1_D(round(x_rspar/dl)) = load1_D(round(x_rspar/dl))+ R_R_D -1/4*n*W_wing_total;
load1_D(round(x_htail/dl)) = load1_D(round(x_htail/dl))- L_tail_D - n*W_htail;
load1_D(round(x_vtail/dl)) = load1_D(round(x_vtail/dl)) - n*W_vtail;
load1_D(round(x_mlg/dl)) = load1_D(round(x_mlg/dl)) - n*W_mlg*2;
load1_D(round(x_nlg/dl)) = load1_D(round(x_nlg/dl)) - n*W_nlg;

% Plot Load Distribution
figure(1)
plot(distance,load1_A,'LineStyle','-');
title("Weight distribution along fuselage");
xlabel("Length (m)");
ylabel("Weight distribution (N)");
% ylim([-100000,100000])
grid on

figure(2)
plot(distance,load1_D,'LineStyle','-');
title("Weight distribution along fuselage");
xlabel("Length (m)");
ylabel("Weight distribution (N)");
% ylim([-10000,10000])
grid on



% LOAD CASE 2: OEI
% create torsion on the fuselage
V_vs = 1.885e5; % shear force on vertical stabliser
l_vs = 6.80*1.3+6/2; % moment arm for torsion
T_vs = 1.885e5*l_vs; % torsion on fuselage

% discretize
load2 = zeros(1000,1);

% find spar reaction
A = [l_F, -l_R;
    1, 1];
B = [-L_tail_D+l2;
    W0+L_tail_D];
x = A^(-1)*B;
R_F_OEI = x(1) ; % reaction force on front spar
R_R_OEI = x(2); % reaction force on rear spar

% populate values
% uniform load due to fuselage empty weight, luggage, passenger
load2(:) = load2(:) -(q_fuselage+q_pass+q_lugg)*dl; 
% load due to fuel tank
load2(round(23.9/dl):round(32.6/dl)) = load1_A(round(23.9/dl):round(32.6/dl)) - q_fuel*dl;
load2(round(x_fspar/dl)) = load2(round(x_fspar/dl))+ R_F_OEI - 3/4*W_wing_total;
load2(round(x_rspar/dl)) = load2(round(x_rspar/dl))+ R_R_OEI - 1/4*W_wing_total;
load2(round(x_htail/dl)) = load2(round(x_htail/dl))- L_tail_D - W_htail;
load2(round(x_vtail/dl)) = load2(round(x_vtail/dl)) - W_vtail;
load2(round(x_mlg/dl)) = load2(round(x_mlg/dl)) - W_mlg*2;
load2(round(x_nlg/dl)) = load2(round(x_nlg/dl)) - W_nlg;

figure(3)
plot(distance,load2,'LineStyle','-');
title("Weight distribution along fuselage");
xlabel("Length (m)");
ylabel("Weight distribution (N)");
% ylim([-10000,10000])
grid on



% LOAD CASE 3: Landing with main gears only
% discretize
load3 = zeros(1000,1);

% load due to landing gear
n = 3; % landing load factor
gear_load = 6.2349e6; % load on one landing gear
load3(round(x_mlg/dl)) = load3(round(x_mlg/dl)) + gear_load*2;
load3(round(x_mlg/dl)) = load3(round(x_mlg/dl)) - n*W_mlg*2;
load3(round(x_nlg/dl)) = load3(round(x_nlg/dl)) - n*W_nlg;


%%% 1. Reaction Force from Tail (neglect, assume tail produce no lift)

%%%%%% i. Lift (neglect, assume wing produce no lift)
%%%%%% ii. Weight
load3(round(x_htail/dl)) = load3(round(x_htail/dl)) - n*W_htail;
load3(round(x_vtail/dl)) = load3(round(x_vtail/dl)) - n*W_vtail;


%%% 2. Reaction Force from Wing 

%%%%%% i. Lift (neglect, assume wing produce no lift)


%%%%%% ii. Weight
load3(round(x_fspar/dl)) = load3(round(x_fspar/dl)) -3/4*n*W_wing_total;
load3(round(x_rspar/dl)) = load3(round(x_rspar/dl)) -1/4*n*W_wing_total;

%%%%%% iii. Thrust (neglect)


%%% 3. Fuselage Weight

% uniform load due to fuselage empty weight, luggage, passenger
load3(:) = load3(:) -n*(q_fuselage+q_pass+q_lugg)*dl; 
% load due to fuel tank
load3(round(23.9/dl):round(32.6/dl)) = load3(round(23.9/dl):round(32.6/dl)) - n*q_fuel*dl;


figure(4)
plot(distance,load3,'LineStyle','-');
title("Weight distribution along fuselage");
xlabel("Length (m)");
ylabel("Weight distribution (N)");
% ylim([-100000,100000])
grid on

figure(5)
plot(distance,load1_D,'LineStyle',':','LineWidth',1.5);
hold on
plot(distance,load1_A,'LineStyle',"-.",'LineWidth',1.5)
plot(distance,load3,'LineStyle',"--",'LineWidth',1.5)
title("Weight distribution along fuselage");
xlabel("Length (m)");
ylabel("Weight distribution (N)");
legend("D","A","Landing")
% ylim([-100000,100000])
grid on

%%
% calculate shear
shear = zeros(n_discrete);
for i = 2:n_discrete
    shear(i) = shear(i-1) - load1_A(i);
end

% Plot Shear Distribution
figure(3)
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
figure(4)
plot(distance,moment);
title("Moment distribution along fuselage");
xlabel("Length (m)");
ylabel("Moment distribution (N)");
grid on






