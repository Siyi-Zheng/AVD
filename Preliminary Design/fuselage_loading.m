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
x_fspar = 33.4; % location of wing front spar
x_rspar = 40.66; % location of wing rear spar
x_ffuelt = 28.6; % location of front fuel tank
x_rfuelt = 37.3;  % location of rear fuel tank

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
l2 = abs(x_htail-xw_AC); % distance of tail AC (assume AC coincide with CG) from aircraft AC
[~,a,~,rho,~] = atmosisa(0);  
M_A = 0.5*rho*(V_A)^2*S*MAC*C_M; % pitching moment at VA
M_D = 0.5*rho*(V_D)^2*S*MAC*C_M; % pitching moment at VD
L_tail_A = (n*W0*l1-M_A)/l2; % lift produce by tail at VA (downwards positive)
L_tail_D = (n*W0*l1-M_D)/l2; % lift produce by tail at VD (downwards positive)

%%% 2. Reaction Force from Wing

%%%%%% i. Lift
L_wing = W0*n; 
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
q_pass = 45*47.88*1.5; % load per unit floor area
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
load1_A(1:round(x_fspar/dl))  = load1_A(1:round(x_fspar/dl)) - 1.18*(q_fuselage+q_pass+q_lugg)*dl;
load1_A(round(x_fspar/dl):end)  = load1_A(round(x_fspar/dl):end) + 0.70*(q_fuselage+q_pass+q_lugg)*dl;
% load due to fuel tanl
load1_A(round(x_ffuelt/dl):round(x_rfuelt/dl)) = load1_A(round(x_ffuelt/dl):round(x_rfuelt/dl)) - q_fuel*n*dl;
load1_A(round(x_fspar/dl)) = load1_A(round(x_fspar/dl)) + R_F_A - 3/4*n*W_wing_total;
load1_A(round(x_rspar/dl)) = load1_A(round(x_rspar/dl)) + R_R_A - 1/4*n*W_wing_total;
load1_A(round(x_htail/dl)) = load1_A(round(x_htail/dl)) - L_tail_A - n*W_htail;
load1_A(round(x_vtail/dl)) = load1_A(round(x_vtail/dl)) - n*W_vtail;
load1_A(round(x_mlg/dl)) = load1_A(round(x_mlg/dl)) - n*W_mlg*2;
load1_A(round(x_nlg/dl)) = load1_A(round(x_nlg/dl)) - n*W_nlg;

% populate load values (VD)
% uniform load due to fuselage empty weight, luggage, passenger
load1_D(:) = load1_D(:) -(q_fuselage+q_pass+q_lugg)*n*dl; 
load1_D(1:round(x_fspar/dl))  = load1_D(1:round(x_fspar/dl)) - 1.62*(q_fuselage+q_pass+q_lugg)*dl;
load1_D(round(x_fspar/dl):end)  = load1_D(round(x_fspar/dl):end) + 1.02*(q_fuselage+q_pass+q_lugg)*dl;
% load due to fuel tanl
load1_D(round(x_ffuelt/dl):round(x_rfuelt/dl)) = load1_D(round(x_ffuelt/dl):round(x_rfuelt/dl)) - q_fuel*n*dl;
load1_D(round(x_fspar/dl)) = load1_D(round(x_fspar/dl))+ R_F_D - 3/4*n*W_wing_total;
load1_D(round(x_rspar/dl)) = load1_D(round(x_rspar/dl))+ R_R_D - 1/4*n*W_wing_total;
load1_D(round(x_htail/dl)) = load1_D(round(x_htail/dl))- L_tail_D - n*W_htail;
load1_D(round(x_vtail/dl)) = load1_D(round(x_vtail/dl)) - n*W_vtail;
load1_D(round(x_mlg/dl)) = load1_D(round(x_mlg/dl)) - n*W_mlg*2;
load1_D(round(x_nlg/dl)) = load1_D(round(x_nlg/dl)) - n*W_nlg;

% % Plot Load Distribution
% figure(1)
% plot(distance,load1_A,'LineStyle','-');
% title("Weight distribution along fuselage");
% xlabel("Length (m)");
% ylabel("Weight distribution (N)");
% % ylim([-100000,100000])
% grid on
% 
% figure(2)
% plot(distance,load1_D,'LineStyle','-');
% title("Weight distribution along fuselage");
% xlabel("Length (m)");
% ylabel("Weight distribution (N)");
% % ylim([-10000,10000])
% grid on



% LOAD CASE 2: OEI 
% equivalent as cruise at 1g + force on vertical stabliser

%%% 1. Torsion on the fuselage
V_vs = 1.885e5; % shear force on vertical stabliser
l_vs = 6.80*1.3+6/2; % moment arm for torsion
T_vs = 1.885e5*l_vs; % torsion on fuselage

% discretize
load2_v = zeros(1000,1);
load2_h = zeros(1000,1);
load2_h(round(x_vtail/dl)) = load2_h(round(x_vtail/dl)) + T_vs;

%%% 2. Lift on Vertical Tailplane
V_OEI = V_A;
M_OEI = 0.5*rho*(V_OEI)^2*S*MAC*C_M; % pitching moment at VD
L_tail_OEI = (W0*l1-M_OEI)/l2; % lift produce by tail at VA (downwards positive)

% find spar reaction
A = [l_F, -l_R;
    1, 1];
B = [-L_tail_OEI+l2;
    W0+L_tail_OEI];
x = A^(-1)*B;
R_F_OEI = x(1) ; % reaction force on front spar
R_R_OEI = x(2); % reaction force on rear spar

% populate values
% uniform load due to fuselage empty weight, luggage, passenger
load2_v(:) = load2_v(:) - 1.1*(q_fuselage+q_pass+q_lugg)*dl; 
load2_v(1:round(x_fspar/dl))  = load2_v(1:round(x_fspar/dl)) - 0.4*(q_fuselage+q_pass+q_lugg)*dl;
load2_v(round(x_fspar/dl):end)  = load2_v(round(x_fspar/dl):end) + 0.4*(q_fuselage+q_pass+q_lugg)*dl;
% load due to fuel tank
load2_v(round(x_ffuelt/dl):round(x_rfuelt/dl)) = load2_v(round(x_ffuelt/dl):round(x_rfuelt/dl)) - 2*q_fuel*dl;
load2_v(round(x_fspar/dl)) = load2_v(round(x_fspar/dl))+ 0.93*R_F_OEI - 3/4*W_wing_total;
load2_v(round(x_rspar/dl)) = load2_v(round(x_rspar/dl))+ R_R_OEI; % - 1/4*W_wing_total;
load2_v(round(x_htail/dl)) = load2_v(round(x_htail/dl))- 0.5*L_tail_D - W_htail;
load2_v(round(x_vtail/dl)) = load2_v(round(x_vtail/dl)) - W_vtail;
load2_v(round(x_mlg/dl)) = load2_v(round(x_mlg/dl)) - W_mlg*2;
load2_v(round(x_nlg/dl)) = load2_v(round(x_nlg/dl)) - W_nlg;

% figure(3)
% yyaxis left
% plot(distance,load2_v,'LineStyle','-','LineWidth',1.5);
% ylabel("Vertical Load Distribution (N)");
% xlabel("Length (m)");
% xlim([0,l_plane])
% ylim([-5e5,20e5])
% yyaxis right
% plot(distance,load2_h,'LineStyle','-.','LineWidth',1.5)
% % title("Weight distribution along fuselage");
% ylabel("Horizontal Load Distribution (N)");
% ylim([-5e6/8,2.5e6])
% legend("Vertical Load","Horizontal Load",Location="northwest")
% grid on



% LOAD CASE 3: Landing with main gears only
% discretize
load3 = zeros(n_discrete,1);

% load due to landing gear
n = 3; % landing load factor
gear_load = 6.2349e6; % load on one landing gear
load3(round(x_mlg/dl)) = load3(round(x_mlg/dl)) + gear_load*4 - n*W_mlg*2;
load3(round(x_nlg/dl)) = load3(round(x_nlg/dl)) - n*W_nlg;


%%% 1. Reaction Force from Tail (neglect, assume tail produce no lift)

%%%%%% i. Lift (neglect, assume wing produce no lift)
%%%%%% ii. Weight
load3(round(x_htail/dl)) = load3(round(x_htail/dl)) - 5*n*W_htail;
load3(round(x_vtail/dl)) = load3(round(x_vtail/dl)) - 5*n*W_vtail;


%%% 2. Reaction Force from Wing 

%%%%%% i. Lift (neglect, assume wing produce no lift)


%%%%%% ii. Weight
load3(round(x_fspar/dl)) = load3(round(x_fspar/dl)) -1.58*3/4*n*W_wing_total;
load3(round(x_rspar/dl)) = load3(round(x_rspar/dl)) -1.02*1/4*n*W_wing_total;

%%%%%% iii. Thrust (neglect)


%%% 3. Fuselage Weight

% uniform load due to fuselage empty weight, luggage, passenger
load3(:) = load3(:) -3.08*n*(q_fuselage+q_pass+q_lugg)*dl; 
% load3(1:round(x_fspar/dl))  = load3(1:round(x_fspar/dl)) - (q_fuselage+q_pass+q_lugg)*dl;
% load3(round(x_fspar/dl):end)  = load3(round(x_fspar/dl):end) + 0.4*(q_fuselage+q_pass+q_lugg)*dl;
% load due to fuel tank
load3(round(x_ffuelt/dl):round(x_rfuelt/dl)) = load3(round(x_ffuelt/dl):round(x_rfuelt/dl)) - n*q_fuel*dl;


% figure(4)
% plot(distance,load3,'LineStyle','-');
% title("Weight distribution along fuselage");
% xlabel("Length (m)");
% ylabel("Weight distribution (N)");
% % ylim([-100000,100000])
% grid on

figure(1)
yyaxis left
plot(distance,load1_D,'LineStyle','-','LineWidth',1.75,'Color',"#0072BD");
hold on
plot(distance,load1_A,'LineStyle',"-.",'LineWidth',1.5,'Color',"r")
plot(distance,load2_v,'LineStyle',"--",'LineWidth',1.5,'Color',"g")
plot(distance,load3,'LineStyle',":",'LineWidth',2,'Color',"#EDB120")
xlabel("Length (m)","FontSize",14);
ylabel("Vertical Load Distribution (N/m)","FontSize",14);
ylim([-30000,30000])
xlim([0,l_plane])
yyaxis right
plot(distance,load2_h,'LineStyle','-','LineWidth',1.5,'Color','m')
% title("Weight distribution along fuselage");
ylabel("Horizontal Load Distribution (N/m)","FontSize",14);
ylim([-2.5e6,2.5e6])

legend("Case 1a","Case 1b","Case 2", "Case 3", "Case 2 - Horizontal","Location","northwest","FontSize",12)
ax = gca;
ax.YAxis(1).Color = '#0072BD';
ax.YAxis(2).Color = 'm';
grid on
hold off

%%
% calculate shear
shear1_D = zeros(n_discrete,1);
for i = 2:n_discrete
    shear1_D(i) = shear1_D(i-1) - dl*load1_D(i);
end

% calculate moment
moment1_D = zeros(n_discrete,1);
for i = 2:n_discrete
    moment1_D(i) = moment1_D(i-1) + dl*shear1_D(i);
end



shear1_A = zeros(n_discrete,1);
for i = 2:n_discrete
    shear1_A(i) = shear1_A(i-1) - dl*load1_A(i);
end

% calculate moment
moment1_A = zeros(n_discrete,1);
for i = 2:n_discrete
    moment1_A(i) = moment1_A(i-1) + dl*shear1_A(i);
end



shear2_v = zeros(n_discrete,1);
for i = 2:n_discrete
    shear2_v(i) = shear2_v(i-1) - dl*load2_v(i);
end
% calculate moment
moment2_v = zeros(n_discrete,1);
for i = 2:n_discrete
    moment2_v(i) = moment2_v(i-1) + dl*shear2_v(i);
end



shear2_h = zeros(n_discrete,1);
for i = 2:n_discrete
    shear2_h(i) = shear2_h(i-1) - dl*load2_h(i);
end
% calculate moment
moment2_h = zeros(n_discrete,1);
for i = 2:n_discrete
    moment2_h(i) = moment2_h(i-1) + dl*shear2_h(i);
end



shear3 = zeros(n_discrete,1);
for i = 2:n_discrete
    shear3(i) = shear3(i-1) - dl*load3(i);
end
% calculate moment
moment3 = zeros(n_discrete,1);
for i = 2:n_discrete
    moment3(i) = moment3(i-1) + dl*shear3(i);
end


%%
% Plot Shear Distribution
figure(2)
yyaxis left
plot(distance,shear1_D);
xlabel("Length (m)");
ylabel("Shear distribution (N/m)");
yyaxis right
plot(distance,moment1_D);
ylabel("Moment distribution (Nm)");
title("CASE-1D")
grid on

% Plot Shear Distribution
figure(3)
yyaxis left
plot(distance,shear1_A);
xlabel("Length (m)");
ylabel("Shear distribution (N/m)");
yyaxis right
plot(distance,moment1_A);
ylabel("Moment distribution (Nm)");
title("CASE-1A");
grid on

% Plot Shear Distribution
figure(4)
yyaxis left
plot(distance,shear2_v);
xlabel("Length (m)");
ylabel("Shear distribution (N/m)");
yyaxis right
plot(distance,moment2_v);
ylabel("Moment distribution (Nm)");
title("CASE-2V");
%ylim([0 15000000])
grid on

% Plot Shear Distribution
figure(5)
yyaxis left
plot(distance,shear2_h);
xlabel("Length (m)");
ylabel("Shear distribution (N/m)");
%ylim([0 150000000])
yyaxis right
plot(distance,moment2_h);
ylabel("Moment distribution (Nm)");
title("CASE-2H");
grid on

% Plot Shear Distribution
figure(6)
yyaxis left
plot(distance,shear3);
xlabel("Length (m)");
ylabel("Shear distribution (N/m)");
yyaxis right
plot(distance,moment3);
ylabel("Moment distribution (Nm)");
title("CASE-3");
grid on

%%
figure(7)
plot(distance,shear1_D,LineWidth=1.5);
hold on 
plot(distance,shear1_A,LineWidth=1.5);
plot(distance,shear2_v,LineWidth=1.5);
plot(distance,shear3,LineWidth=1.5);
xlabel("Length (m)",FontSize=14);
ylabel("Shear distribution (N/m)",FontSize=14);
legend("Case 1a", "Case 1b","Case 2","Case 3",FontSize=14)
ylim([-1e6 1.2e6])
xlim([0 l_plane])
grid on

figure(8)
plot(distance,moment1_D,LineWidth=1.5);
hold on 
plot(distance,moment1_A,LineWidth=1.5);
plot(distance,moment2_v,LineWidth=1.5);
plot(distance,moment3,LineWidth=1.5);
xlabel("Length (m)",FontSize=14);
ylabel("Moment distribution (Nm)",FontSize=14);
legend("Case 1a", "Case 1b","Case 2","Case 3",FontSize=14)
ylim([-2e6 18e6])
xlim([0 l_plane])
grid on






