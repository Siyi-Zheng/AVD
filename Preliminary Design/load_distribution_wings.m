clc
clear
close all

%LOADS DISTRIBUTION

%get chord distribution
%leading edge sweep angle
b = 32.45 * tand(26.6);
a = b - 3.47/4;

aa = (a + 3.47)/32.45;

%CMo for airfoil
Cmoairf = -0.131;

%sweep_angle = atand(aa);
sweep_angle = atand(aa);

%chord distribution for first section
outer_chord = 13.89 - (9.75 * tand(sweep_angle));

gradient1 = (outer_chord - 13.89)/9.75;

%for outer section
gradient2 = (3.47 - outer_chord)/(32.45 - 9.75);

span1 = 3.17:0.01:9.75;
span2 = 9.75:0.01:32.5;

c2 = outer_chord - gradient2 * 9.75;

chord1 = gradient1 * span1 + 13.89;
chord2 = gradient2*span2+c2;

chord = [chord1,chord2];



%constants
W_S = 7840;         %wing loading in N/m^2
S = 482;            %wing area in m^2
rho = 1.227;        %sea level density
W = 353385;         %mtow in kg
g = 9.81;           %gravitational constant

%masses and wing properties
m_dry = 27100/2;        %dry mass of the wing in kg
weight_wing = m_dry*g;  %weight of the wing in N
V_wing = 168.6331;      %volume of the wing in m^3 CHANGE THIS WITH THE REAL VALUE  

m_fuel = 176091;        %fuel mass in kg
rho_fuel = 804;         %fuel density in kg/m^3
length_ft = 26.83;      %length of fuel tank 

m_engine = 6147.1;      %mass of engine
m_mg = 15100;           %mass of main landing gear


%aerofoil properties
cg = 0.4141;                            %centre of gravity
arm_length = 0.1641;                    %distance of cg from aerodynamic centre
frontSpar = 0.12;                       %location of front spar
rearSpar = 0.7;                         %location of shear spar
flexAxis = (frontSpar+rearSpar)/2;      %location of flexural axis (shear centre)


%wing discretisation
span = [span1,span2];      %29.33 is the length of the wing that is outside of the fuselage

%area of the aerofoil at each section of the wing
area_aerofoil = 0.0939*chord.^2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CASE 1 - symmetric flight at ultimate load factor evaluated at manoeuvre
%(speed = V_A)
v_a = 171.25;                  %Manouvre speed
n = 3.75;                        %corresponding load factor from V/N diagram
L = n*W*g;                      %overall lift
C_l = (2*n*W_S)/(rho*v_a^2);    %lift coeff. This is computed to be 1.8749

alpha = 15;     %AoA in degrees taken from polars
C_d = 0.022;     % drag coeff taken from polars

D = 0.5*rho*v_a^2*S*C_d;    %overall drag

%THIS IS ONLY FOR CASE 3, COMMENT OUT OTHERWISE
% L = 0;
% D = 0;
% C_l = 0;
% C_d = 0;

%aerodynamics (assuming an elliptic lift distribution)
L0 = L/(32.5*pi*0.5 - 6.34);

%load on each section of the wing due to the wing's weight
dy = 0.01;                              %spacing between each section
dV = area_aerofoil * dy;                %volume of each section
dW_wing = (weight_wing/V_wing) * dV;    %weight of the wing on each section

%load due to 1 engine in N
l_engine = m_engine*g*n;  

%load due to 1 main landing gear in N
l_mg = m_mg*g*n/4;             

%load due to fuel in the wings in N
V_fuel = 219.75;                %volume of fuel in m^3
dmass_fuel = rho_fuel * dV;     %mass of fuel in each section
dmass_fuel(2685:end) = 0;       %no fuel after y = 30
dW_fuel = dmass_fuel*g;         %weight of fuel in each section

%add all of the weights to get the total weight experienced at each section
total_weight = dW_wing + dW_fuel;                       %add both the loads due to the wing and the fuel
total_weight(623:698) = total_weight(623:698) + l_engine/75;       %add weight of engine
total_weight(1761:1836) = total_weight(1761:1836) + l_engine/75;     %add weight of engine
total_weight(249:297) = total_weight(273) + l_mg/48;           %add weight of main landing gear
total_weight = total_weight*100;
%total_weight(249:297) = total_weight(249:297) - 3.4278e6/48;       %THIS IS ONLY FOR CASE 3 COMMENT OUT OTHERWISE

figure(1)
plot(span,total_weight,"b",LineWidth=1.5);
%title("Weight distribution along the span due to inertia loads");
xlabel("Span (m)");
ylabel("Weight distribution (N)");
grid on


%CALCULATING SHEAR
dS = dy*chord;                      %planform area distribution
dD = 0.5*rho*v_a^2*C_d*dS;          %drag distribution
dL = L0*sqrt(1-(span/32.5).^2);     %elliptic lift distribution

figure(2)
plot(span,dL);
title("Lift distribution");
xlabel("span (m)");
ylabel("Lift (N)");
grid on

for i = 1:length(span)
    temp = dD(i:length(span));
    shear_x(i) = sum(temp*dy);
end

shear_y = 0;

for i = 1:length(span)
    temp1 = dL(i:length(span));
    temp2 = total_weight(i:length(span));
    shear_z(i) = sum((temp1-temp2)*dy);
end

shear_total = sqrt(shear_x.^2 + shear_z.^2);

figure(3)
plot(span,shear_total);
title("shear distribution");
xlabel("span (m)");
ylabel("shear (N/m)");
grid on

%CALCULATING BENDING MOMENT
for i = 1:length(span)
    temp = shear_z(i:length(span));
    M_x(i) = sum(temp);
end

for i = 1:length(span)
    temp = shear_x(i:length(span));
    M_z(i) = sum(temp);
end

figure(4)
plot(span,M_x);
title("x bending moment distribution");
xlabel("span (m)");
ylabel("bending moment (Nm)");
grid on

figure(5)
plot(span,M_z);
title("z bending moment distribution");
xlabel("span (m)");
ylabel("bending moment (Nm)");
grid on


%CALCULATING TORQUE
M_torque = zeros(length(span));
for i=1:length(span)
    x_f(i) = 0.25 * chord(i);
    z_f = 0.07; % depends how we define the axis
    x_cg= cg; 
    dM(i) = -dL(i)* x_f(i) + dD(i) * z_f + total_weight(i) *x_cg + Cmoairf * 3.75;
    M_torque(i)= sum(dM) ;
end

figure(6)
plot(span,shear_z);
title("shear_z distribution");
xlabel("span (m)");
ylabel("shear (N/m)");
grid on

figure(7)
plot(span,shear_x);
title("shear_x distribution");
xlabel("span (m)");
ylabel("shear (N/m)");
grid on



figure(8)
plot(span,dL-total_weight,"b",LineWidth=1.5);
hold on
plot(span,dL,"r");
plot(span,-total_weight);
legend("Load distribution","Lift distribution","Weight distribution");
%title("Weight distribution along the span due to inertia loads");
xlabel("Span (m)");
ylabel("Weight distribution (N/m)");
grid on

figure(9)
plot(span, M_torque)














