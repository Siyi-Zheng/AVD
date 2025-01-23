clear
clc

% Aircraft Parameters
l_plane = 77.8; % length of aircraft
W0 = 354000*9.81; % aircraft weight
W_wing = 27100*9.81; % wing weight
W_engine = 6200*9.81; % engine weight
W_fuel = 1000*9.81; % fuel weight
W_fuselage = 31300*9.81; % fuselage weight
W_tail = (1970+1900)*9.81; % tailplane weight
l_wing = 34.8;
l_tail = 72.1;

% Force on Fuselage

%%% 1. Reaction Force from Wing

%%%%%% i. Lift
n = 2.5; % load factor
L_wing = W0; 

%%%%%% ii. Weight
W_wing_total = W_wing+W_engine+W_fuel;
%%%%%% iii. Thrust


%%% 2. Reaction Force from Tail
[~,a,~,rho,~,~] = atmosisa(10000);
C_M = 0.73;
V_cruise = a*0.83;
l = 72.1-34.8; % distance of tail from aircraft CG
S = 482;
c = 7.41;
L_tail = (0.5*rho*(V_cruise)^2*S*c*C_M)/l; % lift produce by tail


%%% 3. Fuselage Weight
q_fuselage = W_fuselage/l_plane; % load per unit length 

% discretize
n_discrete = 1000;
dl = l_plane/n_discrete;
distance = linspace(0,l_plane,n_discrete);
load = zeros(1000);

% populate load values
load(:) = -q_fuselage*dl;
load(round(l_wing/dl)) = load(round(l_wing/dl))+ L_wing-W_wing_total;
load(round(l_tail/dl)) = load(round(l_tail/dl))+ L_tail-W_tail;

% Plot Load Distribution
figure(1)
plot(distance,load);
title("Weight distribution along fuselage");
xlabel("Length (m)");
ylabel("Weight distribution (N)");
grid on

% calculate shear
shear = zeros(n_discrete);
for i = 2:n_discrete
    shear(i) = shear(i-1) - load(i);
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






