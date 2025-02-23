clear
clc



%%
%selected material properites (probably going to change this!!)
sigma_y = 431*10^6;
E = 73850000000;
G = 28700000000;
density = 2765;
v = 0.3365;

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




%%
% calculate shear
shear1_D = zeros(n_discrete,1);
for i = 2:n_discrete
    shear1_D(i) = shear1_D(i-1) - load1_D(i);
end

% calculate moment
moment1_D = zeros(n_discrete,1);
for i = 2:n_discrete
    moment1_D(i) = moment1_D(i-1) + shear1_D(i);
end

shear1_A = zeros(n_discrete,1);
for i = 2:n_discrete
    shear1_A(i) = shear1_A(i-1) - load1_A(i);
end

% calculate moment
moment1_A = zeros(n_discrete,1);
for i = 2:n_discrete
    moment1_A(i) = moment1_A(i-1) + shear1_A(i);
end

shear2_v = zeros(n_discrete,1);
for i = 2:n_discrete
    shear2_v(i) = shear2_v(i-1) - load2_v(i);
end
% calculate moment
moment2_v = zeros(n_discrete,1);
for i = 2:n_discrete
    moment2_v(i) = moment2_v(i-1) + shear2_v(i);
end

shear2_h = zeros(n_discrete,1);
for i = 2:n_discrete
    shear2_h(i) = shear2_h(i-1) - load2_h(i);
end
% calculate moment
moment2_h = zeros(n_discrete,1);
for i = 2:n_discrete
    moment2_h(i) = moment2_h(i-1) + shear2_h(i);
end

shear3 = zeros(n_discrete,1);
for i = 2:n_discrete
    shear3(i) = shear3(i-1) - load3(i);
end
% calculate moment
moment3 = zeros(n_discrete,1);
for i = 2:n_discrete
    moment3(i) = moment3(i-1) + shear3(i);
end

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
max_stress(:) = sigma_y;

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

skin_thickness = max3 / sigma_y;

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

hoop_stress = (P * D) / (2 * skin_thickness);

long_stress = (P * D) / (4 * skin_thickness);

D_th = D/(2*skin_thickness);
D_tl = D/(4*skin_thickness);

sigma_allow = sigma_y;

stress_ratio = sigma_allow/P;

stress_h = D_th * P;
stress_l = D_tl * P;

%check stress is below maximum allowable stress
assert(stress_h < (sigma_y/3), 'Error: stress is too high');
assert(stress_l < (sigma_y/3), 'Error: stress is too high');

%getting hemispherical ends
thickness_ratio = (2-v)/(1-v);

hemispherical_thickness = skin_thickness / thickness_ratio;
%%
%ADD STRINGER CODE HERE
%constants
moment_x = 2.6393e+09; % bending moment about x (Nm)
moment_y = 1000; % bending moment about y (Nm)
d_fuslg = 6.3; % fuselage diameter (m)
circum_fuslg = pi*d_fuslg; % fuselage circumference (m)
skin_thickness;

%inputs
no_stringer = 80; % number of stringers
dl_stringer = circum_fuslg/no_stringer; % stringer spacing
fprintf("Stringer Spacing : %.3g inches \n",dl_stringer*39.3701)
A_stringer = 60e-6; % one stringer area (m^2)
dA_boom = ones(1,no_stringer) * 15*skin_thickness^2; % skin collaborative area
A_boom =  A_stringer + dA_boom; % boom area (m^2)
fprintf("Boom Area (0) : %.3g \n",A_boom(1));


%calculating boom distribution
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


dA_boom(1) = skin_thickness/6 * (2 + sig(no_stringer)/sig(1)) + skin_thickness/6 * (2 + sig(2)/sig(1));
for i = 2:no_stringer-1
    dA_boom(i) = skin_thickness/6 * (2 + sig(i-1)/sig(i)) + skin_thickness/6 * (2 + sig(i+1)/sig(i));
end
dA_boom(no_stringer) = skin_thickness/6 * (2 + sig(no_stringer-1)/sig(no_stringer)) + skin_thickness/6 * (2 + sig(1)/sig(no_stringer));
A_boom =  A_stringer + dA_boom;
fprintf("Boom Area (1) : %.3g \n",A_boom(2));

for i = 1:5
    Ix = sum(A_boom.*(y.*y)); % second moment of area about x
    Iy = sum(A_boom.*(x.*x)); % second moment of area about y
    sigx = moment_x*y/Ix; % direct stress due to x bending moment
    sigy = moment_y*x/Iy; % direct stress due to y bending moment
    sig = sigx + sigy; % direct stress due to x and y bending moment
    dA_boom(1) = skin_thickness/6 * (2 + sig(no_stringer)/sig(1)) + skin_thickness/6 * (2 + sig(2)/sig(1));
    for j = 2:no_stringer-1
        dA_boom(j) = skin_thickness/6 * (2 + sig(j-1)/sig(j)) + skin_thickness/6 * (2 + sig(j+1)/sig(j));
    end
    dA_boom(no_stringer) = skin_thickness/6 * (2 + sig(no_stringer-1)/sig(no_stringer)) + skin_thickness/6 * (2 + sig(1)/sig(no_stringer));
    A_boom =  A_stringer + dA_boom;
    fprintf("Boom Area (%d) : %.3g \n",i+1,A_boom(2));
end


sig_max = max(sig);

% plot boom location
figure(1)
plot(x,y,".",MarkerSize=8)

% Plot fuselage cross-section
figure(2);
hold on;
grid on;
plot3(x, y, z, 'LineWidth', 2,'Color','r'); % Circular fuselage cross-section
plot3(x, y, z, '.','MarkerSize', 12,'Color','b') 
% Plot stress vectors at each stringer location
quiver3(x, y, z, u, v, sigx+sigy, 'k', 'LineWidth', 1, 'MaxHeadSize', 0.5);
% Labels and title
xlabel('x');
ylabel('y');
zlabel('Direct stress \sigma_z (MPa)');
title('Direct stress distribution around fuselage');
legend({'Fuselage cross-section', 'String Location','Direct Stress'}, 'Location', 'northeast');
view(3); % Set 3D view
hold off;


% using Z shape stringer
syms t_stringer; % stringer thickness
syms w_flange;  % stringer flange width
syms h_web; % stringer web height
Area = 2*w_flange*t_stringer + (h_web - 2*t_stringer)*t_stringer == A_stringer; % total area of stringer
% need to be the same as A_stringer
I = (t_stringer*h_web^3)/12 + 2*(t_stringer*w_flange*(h_web/2)^2); % second moment of area of stringer
f = solve(Area,h_web);
f = matlabFunction(f, 'Vars', [t_stringer, w_flange]);
fsurf(f,[1e-3 5e-3 5e-2 20e-2])


% check for failure
% material: Al 2024-T861
sig_yield = 431e6; % yield stress
fprintf("Maximum stress: %.2g MPa \nYield stress: %.2g MPa",sig_max/1e6,sig_yield/1e6);

% check for stringer buckling
E = 73.85e9; % Young's Modulus
syms L_eff % effective length
sig_euler = pi^2*E*I/(L_eff^2*A_stringer); % euler buckling stress

% check for skin buckling
sig_skin = 3.6*E*(skin_thickness/dl_stringer)^2;








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

h = 0.01:0.01:0.2;

b = 0.01:0.01:0.6;

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
xlabel('Base Width (b)');
ylabel('Frame Height (h)');
zlabel('Frame Area (Af)');
title('Frame Area vs. Base Width and Height');
set(gca,'XDir','reverse','YDir','reverse')
shading interp; % Smooth shading
s.EdgeColor = "[0,0,0]";
colorbar; % Add a color legend


figure
s = surf(H', B', t); % Use correctly shaped variables
xlabel('Base Width (b)');
ylabel('Frame Height (h)');
zlabel('Frame Thickness (t)');
title('Frame Thickness vs. Base Width and Height');
set(gca,'XDir','reverse','YDir','reverse')
shading interp; % Smooth shading
s.EdgeColor = "[0,0,0]";
colorbar; % Add a color legend





%%
clear P Q R
%heavy frames

%Defining key loading points
%need loads from tailplane, landing gear, and wings



