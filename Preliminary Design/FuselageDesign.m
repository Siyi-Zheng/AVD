clear
clc



%%
%selected material properites (probably going to change this!!)
sigma_y = 378*10^6;
E = 73850000000;
G = 28700000000;
density = 2765;
v = 0.3365;
tau_y = sigma_y/sqrt(3);

%%
%load calculation

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

figure('Position', [150, 150, 800*1.1, 300*1.1])
plot(distance,shear1_D,LineWidth=2,Color='r');
hold on 
plot(distance,shear1_A,LineWidth=2,Color='b',LineStyle='--');
plot(distance,shear2_v,LineWidth=2,Color='k',LineStyle='-.');
plot(distance,shear3,LineWidth=2,Color='g',LineStyle=':');
xlabel("Length (m)",FontSize=14);
ylabel("Shear distribution (N/m)",FontSize=14);
legend("Case 1a", "Case 1b","Case 2","Case 3",FontSize=14)
ylim([-1e6 1.2e6])
xlim([0 l_plane])
grid on

figure('Position', [150, 150, 800*1.1, 300*1.1])
plot(distance,moment1_D,LineWidth=2,Color='r');
hold on 
plot(distance,moment1_A,LineWidth=2,Color='b',LineStyle='--');
plot(distance,moment2_v,LineWidth=2,Color='k',LineStyle='-.');
plot(distance,moment3,LineWidth=2,Color='g',LineStyle=':');
xlabel("Length (m)",FontSize=14);
ylabel("Moment distribution (Nm)",FontSize=14);
legend("Case 1a", "Case 1b","Case 2","Case 3",FontSize=14)
ylim([-2e6 18e6])
xlim([0 l_plane])
grid on

max_shear = max(abs(shear3));
max_bending = max(abs(moment3));
max_torsion = T_vs;

fprintf("Max Shear Stress: %.3g \nMax Bending Moment: %.3g \nMax Torsion: %.3g \n",max_shear,max_bending,max_torsion)

%%
%skin thickness shear flow
%need to break load cases down into 

q1A = shearflowplot(0 , shear1_A , 0);
q1D = shearflowplot(0 , shear1_D , 0);
q2 = shearflowplot(shear2_h , shear2_v , T_vs);
q3  = shearflowplot(0 , shear3 , 0);

%max shear flow is needed
[max1A , I1A] = max(abs(q1A(:)));
[max1D , I1D] = max(abs(q1D(:)));
[max2 , I2]= max(abs(q2(:)));
[max3 , I3] = max(abs(q3(:)));

[row_max1A, col_max1A] = ind2sub(size(q1A), I1A);
[row_max1D, col_max1D] = ind2sub(size(q1D), I1D);
[row_max2, col_max2] = ind2sub(size(q2), I2);
[row_max3, col_max3] = ind2sub(size(q3), I3);

q_max = shearflowplot(0 , 1.08*10^6 , 2.23*10^6);
q_max = max(q_max);

%plotting shear for max case
%copy data

qmax = q3(500,:);


%data for fuselage ring
r = 6.34;

%normalize qmax
qmax(:) = abs(qmax(:)./max(qmax)) + r;

phi = 0:10:360;
phi = phi .* (pi/180);
rho = zeros(length(phi) , 1);
rho(:) = r;

%fuselage ring
fus_ring = zeros(length(phi) , 1);
fus_ring(:) = min(qmax);

%plotting shear flow distribution
figure

polarplot(phi , rho , Color='r' , LineWidth=2 , DisplayName='Fusalge Skin')
hold on
polarplot(phi , qmax , Color='c' , LineWidth=2 , DisplayName='Shear Flow')

ax = gca;
d = ax.ThetaDir;
ax.ThetaDir = 'counterclockwise';
ax.ThetaZeroLocation = 'bottom';
legend(FontSize=14)
hold off

%max thickness needed
t = 0:0.0001:0.01;

%take maximum acceptable yield stress as 1/3 of yield stress
max_stress = zeros(length(t) , 1);
max_stress(:) = sigma_y * 0.58;

%t in mm
t_mm = t*1000;

figure
hold on
set(gca, 'YScale', 'log') % Force log scale on Y-axis
semilogy(t_mm , max_stress/1000000 , 'r--')
semilogy(t_mm , (max3./t/1000000))

ylabel('Stress (Mpa)')
xlabel('Skin Thickess (m)')
hold off

skin_thickness = q_max / (sigma_y * 0.58);

%%
%Pressurization Loads
%taking pressure difference at cruise
Ho = 43000;
Hi = 8000;
Ho = convlength(Ho,'ft','m');
Hi = convlength(Hi,'ft','m');
format long g;
[~,~,Po] = atmosisa(Ho);
[~,~,Pi] = atmosisa(Hi);

P = Pi - Po;

%fuselage is perfect cylinder with diameter 
D = 6.34;

stress_hoop = (P * D) / (2 * 0.0015);

stress_long = (P * D) / (4 * 0.0015);

D_th = D/(2*0.001);
D_tl = D/(4*0.001);

sigma_allow = sigma_y/1.5;

skint_required = (P * D) / (2 * sigma_allow);

stress_ratio = sigma_allow;

stress_h = D_th * P;
stress_l = D_tl * P;

%{
%check stress is below maximum allowable stress
assert(stress_h < (sigma_allow), 'Error: stress is too high');
assert(stress_l < (sigma_allow), 'Error: stress is too high');
%}

%getting hemispherical ends
thickness_ratio = (2-v)/(1-v);

hemispherical_thickness = skint_required / thickness_ratio;

total_skin_t = skint_required + skin_thickness;

%calculate stress using 1.5 mil
stress = P * D / 2 * 0.0015;

%%
% Stringer and Direct Stress

no_stringer = 115; % number of stringers
A_stringer = 52e-6; % one stringer area (m^2)
Lfs = 0.60; % light frame separation
t_skin = 0.0015; % skin thickness (m)

moment_x = max_bending; % bending moment about x (Nm)
moment_y = 10; % bending moment about y (Nm)
d_fuslg = 6.3; % fuselage diameter (m)
circum_fuslg = pi*d_fuslg; % fuselage circumference (m)
dl_stringer = circum_fuslg/no_stringer; % stringer spacing
fprintf("Stringer Spacing : %.3g cm, %.3g inches \n",dl_stringer*100,dl_stringer*39.3701)
dA_boom = ones(1,no_stringer) * 15*t_skin*t_skin; % skin collaborative area
A_boom =  A_stringer + dA_boom; % boom area (m^2)
fprintf("Boom Area : %.3g , Stringer Area: %.3g\n",A_boom(1),A_stringer);


boom_angle = linspace(0,360,no_stringer+1);
boom_angle = boom_angle(1:no_stringer);
x = d_fuslg/2 * cos(deg2rad(boom_angle)); % x location of booms
y = d_fuslg/2 * sin(deg2rad(boom_angle)); % y location of booms
z = zeros(1,no_stringer);
u = zeros(1,no_stringer);
v = zeros(1,no_stringer);

Ix = sum(A_boom.*(y.*y)); % second moment of area about x
Iy = sum(A_boom.*(x.*x)); % second moment of area about y

x_max = max(x); % max x distance from centre
y_max = max(y); % max y distance from centre

sigx = moment_x*y/Ix; % direct stress due to x bending moment
sigy = moment_y*x/Iy; % direct stress due to y bending moment
sig = sigx + sigy; % direct stress due to x and y bending moment


dA_boom(1) = t_skin/6 * (2 + sig(no_stringer)/sig(1)) + t_skin/6 * (2 + sig(2)/sig(1));
for i = 2:no_stringer-1
    dA_boom(i) = t_skin/6 * (2 + sig(i-1)/sig(i)) + t_skin/6 * (2 + sig(i+1)/sig(i));
end
dA_boom(no_stringer) = t_skin/6 * (2 + sig(no_stringer-1)/sig(no_stringer)) + t_skin/6 * (2 + sig(1)/sig(no_stringer));
A_boom =  A_stringer + dA_boom;
fprintf("Boom Area (1) : %.3g \n",A_boom(2));

for i = 1:5
    Ix = sum(A_boom.*(y.*y)); % second moment of area about x
    Iy = sum(A_boom.*(x.*x)); % second moment of area about y
    sigx = moment_x*y/Ix; % direct stress due to x bending moment
    sigy = moment_y*x/Iy; % direct stress due to y bending moment
    sig = sigx + sigy; % direct stress due to x and y bending moment
    dA_boom(1) = t_skin/6 * (2 + sig(no_stringer)/sig(1)) + t_skin/6 * (2 + sig(2)/sig(1));
    for j = 2:no_stringer-1
        dA_boom(j) = t_skin/6 * (2 + sig(j-1)/sig(j)) + t_skin/6 * (2 + sig(j+1)/sig(j));
    end
    dA_boom(no_stringer) = t_skin/6 * (2 + sig(no_stringer-1)/sig(no_stringer)) + t_skin/6 * (2 + sig(1)/sig(no_stringer));
    A_boom =  A_stringer + dA_boom;
    fprintf("Boom Area (%d) : %.3g \n",i+1,A_boom(2));
end

sig_max = max(sig);


% Plot fuselage cross-section
figure();
hold on;
grid on;
plot3(x, y, z, 'LineWidth', 2,'Color','r'); % Circular fuselage cross-section
plot3(x, y, z, '.','MarkerSize', 12,'Color','b') 
% Plot stress vectors at each stringer location
quiver3(x, y, z, u, v, sigx+sigy, 'k', 'LineWidth', 1, 'MaxHeadSize', 0.5);
% Labels and title
xlabel('x',FontSize=14);
ylabel('y',FontSize=14);
zlabel('Direct stress \sigma_z (MPa)',FontSize=14);
% title('Direct stress distribution around fuselage');
legend({'Fuselage cross-section', 'String Location','Direct Stress'}, 'Location', 'northeast',FontSize=14);
view(3); 

% using Z shape stringer
syms t_stringer; % stringer thickness
syms w_flange;  % stringer flange width
syms h_web; % stringer web height
Area = 2*w_flange*t_stringer + (h_web - 2*t_stringer)*t_stringer == A_stringer; % total area of stringer, need to be the same as A_stringer
I = (t_stringer*h_web^3)/12 + 2*(t_stringer*w_flange*(h_web/2)^2); % second moment of area of stringer


% Define the objective function (negative for maximization)
objective = @(x) -((x(1)*x(2)^3)/12 + 2*(x(1)*x(3)*(x(2)/2)^2));

% Define the constraint function (equality constraint)
constraint = @(x) deal([], 2*x(3)*x(1) + (x(2) - 2*x(1))*x(1) - A_stringer);

% Set initial guesses [t_stringer, h_web, w_flange]
x0 = [2e-3, 0.02, 0.01];

% Define bounds [t_stringer, h_web, w_flange]
lb = [1e-3, 0.01, 0.005]; % Lower bounds (A = 18e-6)
ub = [3e-3, 0.04, 0.02]; % Upper bounds (A = 222e-6)

% Run optimization
options = optimoptions('fmincon', 'Algorithm', 'sqp');
[x_opt, I_max] = fmincon(objective, x0, [], [], [], [], lb, ub, constraint, options);
opt_t_stringer = x_opt(1);
opt_h_web = x_opt(2);
opt_w_flange = x_opt(3);

area = @(x) 2*x(3)*x(1) + (x(2) - 2*x(1))*x(1) - A_stringer;
fprintf("Area constraint: %.2g \n",area(x_opt));


% Display results
disp('Optimal t_stringer, h_web, w_flange:');
disp(x_opt);
disp('Maximum I:');
disp(-I_max); % Since we minimized -I, negate it to get max value

D = 6.34; %frame diameter
Cf = 1/16000; % empirical based off work performed by Lockheed Martin and Shanley, J. Aero Science
M_ult = 1.6e7; % max bending stress
E = 73.85e9; % Young's Modulus

%set a fixed IF_xx
If = (Cf * M_ult * D^2) / (E * Lfs);
% Lfs = (Cf * M_ult * D^2) / (E * If);

%C design 
% If =  t * ((h^3/12) + (0.5 * b * h^2)); % subject to constraint
% Af = (2 * b + h) * t; % minimise area(mass)


% Define the objective function (minimizing Af)
objective = @(x) (2 * x(2) + x(1)) * x(3);

% Define the constraint function (equality constraint)
constraint = @(x) deal([], x(3) * ((x(1)^3 / 12) + (0.5 * x(2) * x(1)^2)) - If);

% Set initial guesses [h, b, t]
x0 = [0.15, 0.03, 2e-3];

% Define bounds [h, b, t]
lb = [0.10, 0.02, 1e-3]; % Lower bounds (I = 1.833e-7, 2.97m)
ub = [0.20, 0.08, 4e-3]; % Upper bounds  (I = 9.067e-6, 0.06m)

% Run optimization
options = optimoptions('fmincon', 'Algorithm', 'sqp');
[x_opt, Af_min] = fmincon(objective, x0, [], [], [], [], lb, ub, constraint, options);
I_area = @(x) x(3) * ((x(1)^3 / 12) + (0.5 * x(2) * x(1)^2)) - If;
fprintf("I constraint: %.2g \n",I_area(x_opt));

% Display results
disp('Optimized values:');
disp(['h = ', num2str(x_opt(1))]);
disp(['b = ', num2str(x_opt(2))]);
disp(['t = ', num2str(x_opt(3))]);
disp(['Minimum Af = ', num2str(Af_min)]);


% check for failure
% material: Al 2024-T861
sig_yield = 431e6; % yield stress
fprintf("Maximum stress: %.2g MPa, Yield stress: %.2g MPa \n",sig_max/1e6,sig_yield/1e6);

% check for stringer buckling
% L_eff = [0.8, 1.6, 3.2, 4.8, 6.4]; % effective length 
sig_euler = pi^2*E*(-I_max)/(Lfs^2*A_stringer); % euler buckling stress
% sig_euler = double(subs(sig_euler, [t_stringer, w_flange, h_web, L_eff], [1.5e-3, 10e-3, 9e-3, 0.4]));
fprintf("Stringer buckling stress: %.2g MPa\n" ,sig_euler/1e6);

% check for skin buckling
w = 1.0*t_skin*sqrt(E/sig_max);
if w < dl_stringer*0.8
    sig_skin = 3.6*E*(t_skin/(dl_stringer - w))^2;
else
    sig_skin = 3.6*E*(t_skin/(dl_stringer))^2;
end
fprintf("Skin buckling stress: %.2g MPa \n" ,sig_skin/1e6);


%%

% objective
objective = @(x) mass_function(x(1), x(2), x(3), x(4));  % x = [no_stringer, A_stringer, Lfs, t_skin]

% Constraint function
constraint = @(x) fuslg_constraints(x);  

% Set initial guesses [no_stringer, A_stringer, Lfs, t_skin]
x0 = [100, 100e-6, 0.4, 3e-3];

% Define bounds [no_stringer, A_stringer, Lfs, t_skin]
lb = [50, 200e-6, 0.2, 2.2e-3]; 
ub = [120, 372e-6, 0.5, 4e-3]; 

intcon = 1; 

% Run optimization using genetic algorithm
options = optimoptions('ga', ...
    'Display', 'iter', ...
    'MaxGenerations', 100, ...
    'PopulationSize', 100, ...
    'PlotFcn', @gaplotbestf);  % Built-in function to plot the best objective value


[x_opt, mass_min] = ga(objective, 4, [], [], [], [], lb, ub, constraint, intcon, options);
no_stringer = x_opt(1);
A_stringer = x_opt(2);
Lfs = x_opt(3);
t_skin = x_opt(4);





%%
%light frame design

%need to vary frame spacing, Lfs, frame section shape Sf and Frame section
%dimensions to mnimise light frame mass, mlf

%creating iterative function to vary through h, Lfs and 

%set section shape -- true if section shape is rectanguler
%                  -- false if section is C-shape

shapefactor = false;

%range for Lfs

%Lfs = 0.1:0.1:10;

%range for h
%{
%section shape consideration
% for a fixed I_xx
%comparing areas for both
sectionshape = true;
Lfs = 0.01:0.01:1;
h = 0.15;

for idx = 1:length(Lfs)
[mlf_rec , If_rec , ~ , Af_rec , ~] = lightframes(Lfs,h,b,sectionshape,E,rho);
end
sectionshape = false;
b = 0.02;
for idx = 1:length(Lfs)
[mlf_C , If_C , ~ , Af_C , ~] = lightframes(Lfs,h,b,sectionshape,E,rho);
end

figure
hold on
plot(If_rec , Af_rec , LineWidth=2 , DisplayName='Rectangular Section')
plot(If_C , Af_C , LineWidth=2 , DisplayName='C-Section')
legend;
xlabel('Moment of Inertia of Frame')
ylabel('Required CSA (m^2)')
hold off
%}

h = 0.02:0.01:0.2;

b = 0.02:0.01:0.2;

%nested for loop to iterate for mass

%initialising results
masses = [];
I_xx = [];
nf = [];
Af = [];
t = [];

Lfs = 0.6;

for idx = 1:length(b)
    for idx2 = 1:length(h)
        
        [masses(idx,idx2) , I_xx(idx,idx2) , nf(idx,idx2) , Af(idx,idx2) , t(idx,idx2)] = lightframes(Lfs , h(idx2) , b(idx) ,  shapefactor , E , density);

    end
end



%finding minimum mass for optimal frame seperation and shape
[min_mass , I_min] = min(masses(:));
[row_idx , col_idx] = ind2sub(size(masses), I_min);

[B,H] = meshgrid(b,h);

% Surface plot
figure
s = surf(H', B', Af); % Use correctly shaped variables
ylabel('Base Width (m)');
xlabel('Frame Height (m)');
zlabel('Frame Area (Af)');
%title('Frame Area vs. Base Width and Height');
set(gca,'XDir','reverse','YDir','reverse')
shading interp; % Smooth shading
s.EdgeColor = "k";
colormap(hot(20));
%colorbar; 


figure
s = surf(H', B', t); % Use correctly shaped variables
ylabel('Base Width (m)');
xlabel('Frame Height (m)');
zlabel('Frame Thickness (m)');
%title('Frame Thickness vs. Base Width and Height');
set(gca,'XDir','reverse','YDir','reverse')
shading interp; % Smooth shading
s.EdgeColor = "k";
colormap(hot(20));
%colorbar;




%%

%total mass calcluation
%skin
syms length_fus n_frames

mass_skin = 0.0014 * pi * 6.34 * 77 * density


%stringers 
mass_stringers = no_stringer * A_stringer * 77 * density

%light frames
mass_lf = 0.00024 * density * pi * 6.34 * (77/0.5)

total_fuselage_mass = mass_skin + mass_stringers + mass_lf

%%
clear P Q R
%heavy frames

%Defining key loading points
%need loads from tailplane, landing gear, and wings
%spar loads

%defining critical point loads
Landing_gear = 247*10^5;

Front_Spar = 75.8*10^5;
Rear_spar = 35.7*10^5;

Horizontal_stabiliser_front = 2.27*10^5;
Horizontal_stabiliser_rear = 2.27*10^5;

HS_toruqe = T_vs;



%considering wing spars
Loads = [...
    33.25, 30.5, 22.5, 0;  % Front Wing Spar
    40.57, 14.4, 10.6, 0;  % Rear Wing Spar
    38.86, 67.8, 103.2, 0;  % Landing Gear
    71.00, 1.91, 1.24, 0;  % Horizontal Stabiliser
    71.18, 1.885, 0, 22.3];  % Vertical Stabiliser

%calculating WISE plots for each

Loads = Loads * 10^5;

%front spar wing 
[M_fs1 , N_fs1 , S_fs1 , ~] = heavyframes(Loads(1,2) , Loads(1,3) , Loads(1,4));

[M_fs2 , N_fs2 , S_fs2 , theta] = heavyframes(Loads(1,2) , Loads(1,3) , Loads(1,4));

%correct frame of reference
[M_fs1new , N_fs1new , S_fs1new] = anglechange(M_fs1 , N_fs1 , S_fs1 , 53.6);
[M_fs2new , N_fs2new , S_fs2new] = anglechange(M_fs2 , N_fs2 , S_fs2 , (360-53.6));

%superpose both loads
M_fs = M_fs1new + M_fs2new;
N_fs = N_fs1new + N_fs2new;
S_fs = S_fs1new + S_fs2new;

theta = theta * (180/pi);

%or doing it with both loads assumed to act at the same point
M_fscst = M_fs1 + M_fs2;
N_fscst = N_fs1 + N_fs2;
S_fscst = S_fs1 + S_fs2;


figure
hold on
plot(theta , M_fs , LineWidth=1.5 , DisplayName='Moment')
plot(theta , N_fs , LineWidth=1.5 , DisplayName='Normal Force')
plot(theta , S_fs , LineWidth=1.5 , DisplayName='Shear Force')
legend;
xlabel('Circumferential location (degrees)')
ylabel('Load (NM or N)')
hold off

figure
hold on
plot(theta , M_fscst , LineWidth=1.5 , DisplayName='Moment')
plot(theta , N_fscst , LineWidth=1.5 , DisplayName='Normal Force')
plot(theta , S_fscst , LineWidth=1.5 , DisplayName='Shear Force')
legend;
xlabel('Circumferential location (degrees)')
ylabel('Load (NM or N)')
hold off


%rear spar wing 
[M_rs1 , N_rs1 , S_rs1 , ~] = heavyframes(Loads(2,2) , Loads(2,3) , Loads(2,4));

[M_rs2 , N_rs2 , S_rs2 , ~] = heavyframes(Loads(2,2) , Loads(2,3) , Loads(2,4));

%correct frame of reference
[M_rs1new , N_rs1new , S_rs1new] = anglechange(M_rs1 , N_rs1 , S_rs1 , 53.6);
[M_rs2new , N_rs2new , S_rs2new] = anglechange(M_rs2 , N_rs2 , S_rs2 , (360-53.6));

%superpose both loads
M_rs = M_rs1new + M_rs2new;
N_rs = N_rs1new + N_rs2new;
S_rs = S_rs1new + S_rs2new;

%or doing it with both loads assumed to act at the same point
M_rscst = M_rs1 + M_rs2;
N_rscst = N_rs1 + N_rs2;
S_rscst = S_rs1 + S_rs2;

figure
hold on
plot(theta , M_rs , LineWidth=1.5 , DisplayName='Moment')
plot(theta , N_rs , LineWidth=1.5 , DisplayName='Normal Force')
plot(theta , S_rs , LineWidth=1.5 , DisplayName='Shear Force')
legend;
xlabel('Circumferential location (degrees)')
ylabel('Load (NM or N)')
hold off

figure
hold on
plot(theta , M_rscst , LineWidth=1.5 , DisplayName='Moment')
plot(theta , N_rscst , LineWidth=1.5 , DisplayName='Normal Force')
plot(theta , S_rscst , LineWidth=1.5 , DisplayName='Shear Force')
legend;
xlabel('Circumferential location (degrees)')
ylabel('Load (NM or N)')
hold off


%landing gear 
[M_lg1 , N_lg1 , S_lg1 , ~] = heavyframes(Loads(3,2) , Loads(3,3) , Loads(3,4));

[M_lg2 , N_lg2 , S_lg2 , ~] = heavyframes(Loads(3,2) , Loads(3,3) , Loads(3,4));

%correct frame of reference
[M_lg1new , N_lg1new , S_lg1new] = anglechange(M_lg1 , N_lg1 , S_lg1 , 33.3);
[M_lg2new , N_lg2new , S_lg2new] = anglechange(M_lg2 , N_lg2 , S_lg2 , (360-33.3));

%superpose both loads
M_lg = M_lg1new + M_lg2new;
N_lg = N_lg1new + N_lg2new;
S_lg = S_lg1new + S_lg2new;


%or doing it with both loads assumed to act at the same point
M_lgcst = M_lg1 + M_lg2;
N_lgcst = N_lg1 + N_lg2;
S_lgcst = S_lg1 + S_lg2;


figure
hold on
plot(theta , M_lg , LineWidth=1.5 , DisplayName='Moment')
plot(theta , N_lg , LineWidth=1.5 , DisplayName='Normal Force')
plot(theta , S_lg , LineWidth=1.5 , DisplayName='Shear Force')
legend;
xlabel('Circumferential location (degrees)')
ylabel('Load (NM or N)')
hold off

figure
hold on
plot(theta , M_lgcst, LineWidth=1.5 , DisplayName='Moment')
plot(theta , N_lgcst, LineWidth=1.5 , DisplayName='Normal Force')
plot(theta , S_lgcst, LineWidth=1.5 , DisplayName='Shear Force')
legend;
xlabel('Circumferential location (degrees)')
ylabel('Load (NM or N)')
hold off



%Horizontal Stabiliser 
[M_hs1 , N_hs1 , S_hs1 , ~] = heavyframes(Loads(4,2) , Loads(4,3) , Loads(4,4));

[M_hs2 , N_hs2 , S_hs2 , ~] = heavyframes(Loads(4,2) , Loads(4,3) , Loads(4,4));

%correct frame of reference
[M_hs1new , N_hs1new , S_hs1new] = anglechange(M_hs1 , N_hs1 , S_hs1 , (180-57.1));
[M_hs2new , N_hs2new , S_hs2new] = anglechange(M_hs2 , N_hs2 , S_hs2 , (180+57.1));

%superpose both loads
M_hs = M_hs1new + M_hs2new;
N_hs = N_hs1new + N_hs2new;
S_hs = S_hs1new + S_hs2new;


%or doing it with both loads assumed to act at the same point
M_hscst = M_hs1 + M_hs2;
N_hscst = N_hs1 + N_hs2;
S_hscst = S_hs1 + S_hs2;



figure
hold on
plot(theta , M_hs , LineWidth=1.5 , DisplayName='Moment')
plot(theta , N_hs , LineWidth=1.5 , DisplayName='Normal Force')
plot(theta , S_hs ,  LineWidth=1.5 , DisplayName='Shear Force')
legend;
xlabel('Circumferential location (degrees)')
ylabel('Load (NM or N)')
hold off

figure
hold on
plot(theta , M_hscst ,  LineWidth=1.5 , DisplayName='Moment')
plot(theta , N_hscst ,  LineWidth=1.5 , DisplayName='Normal Force')
plot(theta , S_hscst ,  LineWidth=1.5 , DisplayName='Shear Force')
legend;
xlabel('Circumferential location (degrees)')
ylabel('Load (NM or N)')
hold off




%vertical Stabiliser 
[M_vs1 , N_vs1 , S_vs1 , ~] = heavyframes(Loads(5,2) , Loads(5,3) , Loads(5,4));


%correct frame of reference
[M_vs1new , N_vs1new , S_vs1new] = anglechange(M_vs1 , N_vs1 , S_vs1 , (180));

%superpose both loads
M_vs = M_vs1new ;
N_vs = N_vs1new ;
S_vs = S_vs1new ;


%or doing it with both loads assumed to act at the same point
M_vscst = M_vs1;
N_vscst = N_vs1;
S_vscst = S_vs1;



figure
hold on
plot(theta , M_vs ,  LineWidth=1.5 , DisplayName='Moment')
plot(theta , N_vs ,  LineWidth=1.5 , DisplayName='Normal Force')
plot(theta , S_vs ,  LineWidth=1.5 , DisplayName='Shear Force')
legend;
xlabel('Circumferential location (degrees)')
ylabel('Load (NM or N)')
legend;
hold off

figure
hold on
plot(theta , M_vscst ,  LineWidth=1.5 , DisplayName='Moment')
plot(theta , N_vscst,  LineWidth=1.5 , DisplayName='Normal Force')
plot(theta , S_vscst ,  LineWidth=1.5 , DisplayName='Shear Force')
legend;
xlabel('Circumferential location (degrees)')
ylabel('Load (NM or N)')
hold off

%obtaining maximum loads in each case
Max_M_fs = max(M_fs);
Max_M_fsc = max(M_fscst);
Max_N_fs = max(N_fs);
Max_N_fsc = max(N_fscst);
Max_S_fs = max(S_fs);
Max_S_fsc = max(S_fscst);

Max_M_rs = max(M_rs);
Max_M_rsc = max(M_rscst);
Max_N_rs = max(N_rs);
Max_N_rsc = max(N_rscst);
Max_S_rs = max(S_rs);
Max_S_rsc = max(S_rscst);

Max_M_lg = max(M_lg);
Max_M_lg = max(M_lgcst);
Max_N_lg = max(N_lg);
Max_N_lg = max(N_lgcst);
Max_S_lg = max(S_lg);
Max_S_lg = max(S_lgcst);

Max_M_hs = max(M_hs);
Max_M_hs = max(M_hscst);
Max_N_hs = max(N_hs);
Max_N_hs = max(N_hscst);
Max_S_hs = max(S_hs);
Max_S_hs = max(S_hscst);

Max_M_vs = max(M_vs);
Max_M_vs = max(M_vscst);
Max_N_vs = max(N_vs);
Max_N_vs = max(N_vscst);
Max_S_vs = max(S_vs);
Max_S_vs = max(S_vscst);


%required section areas

%front spar








