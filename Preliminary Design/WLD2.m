clc
clear
close all

%LOADS DISTRIBUTION

%get chord distribution
%leading edge sweep angle
b = 32.45 * tand(26.6);
a = b - 3.47/4;

aa = (a + 3.47)/32.45;

%sweep_angle = atand(aa);
sweep_angle = atand(aa);

%chord distribution for first section
outer_chord = 13.89 - (9.75 * tand(sweep_angle));

gradient1 = (outer_chord - 13.89)/9.75;

%for outer section
gradient2 = (3.47 - outer_chord)/(32.45 - 9.75);

span1 = 2.55:0.01:9.75;
span2 = 9.75:0.01:32.5;

c2 = outer_chord - gradient2 * 9.75;

chord1 = gradient1 * span1 + 13.89;
chord2 = gradient2*span2+c2;

chord = [chord1,chord2];

%CMo for airfoil
Cmoairf = -0.131;

%constants
W_S = 7840;         %wing loading in N/m^2
S = 482;            %wing area in m^2
rho = 1.227;        %sea level density
W = 353385;         %mtow in kg
g = 9.81;           %gravitational constant
T = 321600;         %thrust produced by each engine in N

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


%cases: [case 1a, case 1b, case 3]

v = [130.537; 171.25; 81.21];           %speed in m/s
n = [2.5; 2.5; 3];                      %load factor
L = n*W*g;                              %overall lift
L(3) = 0;                               %assume no aerodynamic forces in case 3
C_l = (2*n*W_S)./(rho*v.^2);            %lift coefficient
C_l(3) = 0;                             %assume no aerodynamic forces in case 3
%C_d = [0.13; 0.022; 0];                 %drag coefficient taken from polars
C_d = [0.043; 0.017; 0];
D = 0.5.*rho.*v.^2.*S.*C_d;             %overall drag
L0 = L./(32.5*pi*0.5 - 6.34);           %assume an elliptical lift distribution

dy = 0.01;                      %length of each section
dV = area_aerofoil * dy;        %volume of each section

%load on each section of the wing due to the wing's weight
dW_wing = (weight_wing/V_wing) * dV;

%load due to 1 engine in N
dW_engine = [m_engine*g;m_engine*g;m_engine*g]; 

%load due to 1 main landing gear in N
dW_mg = [m_mg*g/4;m_mg*g/4;m_mg*g/4];

%load due to fuel in the wings in N
V_fuel = 219.75;                %volume of fuel in m^3
dmass_fuel = rho_fuel * dV;     %mass of fuel in each section
dmass_fuel(2747:end) = 0;       %no fuel after y = 30
dW_fuel = dmass_fuel*g;         %weight of fuel in each section

%add all of the weights to get the total weight experienced at each section
dW_total = zeros(3, length(chord));             %pre-define array

%case 1 - worst case scenario with no fuel weight
dW_total(1,:) = dW_total(1,:) + dW_wing;                             %add weight of wing
dW_total(1,685:760) = dW_total(1,685:760) + dW_engine(1)/75;         %add weight of engine 1
dW_total(1,1823:1898) = dW_total(1,1823:1898) + dW_engine(1)/75;     %add weight of engine 2
dW_total(1,311:359) = dW_total(1,311:359) + dW_mg(1)/48;             %add weight of landing gear
dW_total(1,:) = dW_total(1,:) * 100;

%case 2 - worst case scenario with no fuel weight
dW_total(2,:) = dW_total(2,:) + dW_wing;                             %add weight of wing
dW_total(2,685:760) = dW_total(2,685:760) + dW_engine(2)/75;         %add weight of engine 1
dW_total(2,1823:1898) = dW_total(2,1823:1898) + dW_engine(2)/75;     %add weight of engine 2
dW_total(2,311:359) = dW_total(2,311:359) + dW_mg(2)/48;             %add weight of landing gear
dW_total(2,:) = dW_total(2,:) * 100;

%case 3 - worst case scenario with fuel weight
dW_total(3,:) = dW_total(3,:) + dW_wing + dW_fuel;                   %add weight of wing and fuel
dW_total(3,685:760) = dW_total(3,685:760) + dW_engine(3)/75;         %add weight of engine 1
dW_total(3,1823:1898) = dW_total(3,1823:1898) + dW_engine(3)/75;     %add weight of engine 2
dW_total(3,311:359) = dW_total(3,311:359) + dW_mg(3)/48;             %add weight of landing gear
dW_total(3,:) = dW_total(3,:) * 100;
dW_total(3,311:359) = dW_total(3,311:359) - 6.2349e6/48;             %add force of main landing gear on the wing at touch down

dW_total = dW_total.*n.*1.5;


%loading lifting distribution
data = load("WingLoad.mat");
dL = data.vq;
dL(1 , :) = dL';
dL(2 , :) = dL(1,:);
dL(3 , :) = zeros(1 , 2997);



%CALCULATING SHEAR
dS = dy*chord;                              %planform area distribution
dD = 0.5*rho*v.^2.*C_d.*dS;                 %drag distribution
%dL = L0.*sqrt(1-(span/32.5).^2);            %elliptical lift distribution
dD(:,685:760) = dD(:,685:760) - T/75;       %subtract thrust of engine 1
dD(:,1823:1898) = dD(:,1823:1898) - T/75;   %subtract thrust of engine 2

for i = 1:3  
    for j = 1:length(span)
        temp = dD(i,j:length(span));
        shear_x(i,j) = sum(temp*dy);
    end
end

shear_y = 0;

for i = 1:3
    for j = 1:length(span)
        temp1 = dL(i,j:length(span));
        temp2 = dW_total(i,j:length(span));
        shear_z(i,j) = sum((temp1-temp2)*dy);
    end
end

shear_total = sqrt(shear_x.^2 + shear_z.^2);


%CALCULATING BENDING MOMENT

for i = 1:3
    for j = 1:length(span)
        temp = shear_z(i,j:length(span));
        M_x(i,j) = sum(temp)*dy;
    end
end

for i = 1:3
    for j = 1:length(span)
        temp = shear_x(i,j:length(span));
        M_z(i,j) = sum(temp)*dy;
    end
end


%CALCULATING TORQUE
M_torque = zeros(3, length(span));
landing_gear_force_case3 = zeros(length(span));
landing_gear_force_case3(311:359) = 3.4278e6/48 ; %double check this is correct or not
for i=1:3
    for j = 1: length(span)
        if i==3 % the landing gear will exert a moment
        x_quart(i, j) = 0.25 * chord(j);
        z_f = 0.07; % depends how we define the axis
        x_cg= cg;
        m_0w(i,j) = 0.5* rho * v(i)^2 * chord(j)^2 * Cmoairf; 
        dM(i,j) = dL(i,j)* (flexAxis - x_quart(i,j))+dD(i,j) * 0 + dW_total(i,j) *(flexAxis - x_cg) + landing_gear_force_case3(1,j)*(flexAxis - rearSpar) - m_0w(i,j);
        M_torque(i,j)= sum(dM(i,j)) ;  
        else
    x_quart(i, j) = 0.25 * chord(j);
    z_f = 0.07; % depends how we define the axis
    x_cg= cg;
    m_0w(i,j) = 0.5* rho * v(i)^2 * chord(j)^2 * Cmoairf; 
    dM(i,j) = dL(i,j)* (flexAxis - x_quart(i,j))+dD(i,j) * 0 + dW_total(i,j) *(flexAxis - x_cg) + - m_0w(i,j);
    M_torque(i,j)= sum(dM(i,j)) ;
    end
    end
end


%PLOTSS

figure(1)
plot(span,-dW_total(1,:),"r", LineWidth = 1.5);
hold on
plot(span,-dW_total(2,:),"b--", LineWidth = 1.5);
plot(span,-dW_total(3,:),":", LineWidth = 1.5);
hold off
grid on
xlabel("Semi-span (m)");
ylabel("Weight distribution (N/m)");
legend("Case 1a","Case 1b", "Case 3");

figure(2)
plot(span,dL(1,:),"r", LineWidth = 1.5);
hold on
plot(span,dL(2,:),"b--", LineWidth = 1.5);
plot(span,dL(3,:),":", LineWidth = 1.5);
hold off
title("Lift distribution");
xlabel("Semi-span (m)");
ylabel("Lift (N/m)");
grid on
legend("Case 1a","Case 1b", "Case 3");

figure(3)
plot(span,dL(1,:)-dW_total(1,:),"r", LineWidth = 1.5);
hold on
plot(span,dL(2,:)-dW_total(2,:),"b--", LineWidth = 1.5);
plot(span,dL(3,:)-dW_total(3,:),":", LineWidth = 1.5);
hold off
xlabel("Semi-span (m)");
ylabel("load distribution (N/m)");
grid on
legend("Case 1a","Case 1b", "Case 3");

figure(4)
subplot(1,3,1);
plot(span,shear_x(1,:),"r", LineWidth = 1.5);
hold on
plot(span,shear_x(2,:),"b--", LineWidth = 1.5);
plot(span,shear_x(3,:),":", LineWidth = 1.5);
hold off
grid on
xlabel("Semi-span (m)");
ylabel("Shear distribution in x-direction (N/m^2)");
legend("Case 1a","Case 1b", "Case 3");

subplot(1,3,2);
plot(span,shear_z(1,:),"r", LineWidth = 1.5);
hold on
plot(span,shear_z(2,:),"b--", LineWidth = 1.5);
plot(span,shear_z(3,:),":", LineWidth = 1.5);
hold off
grid on
xlabel("Semi-span (m)");
ylabel("Shear distribution in z-direction (N/m^2)");
legend("Case 1a","Case 1b", "Case 3");

subplot(1,3,3);
plot(span,shear_total(1,:),"r", LineWidth = 1.5);
hold on
plot(span,shear_total(2,:),"b--", LineWidth = 1.5);
plot(span,shear_total(3,:),":", LineWidth = 1.5);
hold off
grid on
xlabel("Semi-span (m)");
ylabel("Total shear distribution (N/m^2)");
legend("Case 1a","Case 1b", "Case 3");

figure(5)
subplot(1,2,1);
plot(span,M_x(1,:),"r", LineWidth = 1.5);
hold on
plot(span,M_x(2,:),"b--", LineWidth = 1.5);
plot(span,M_x(3,:),":", LineWidth = 1.5);
hold off
grid on
xlabel("Semi-span (m)");
ylabel("Bending moment distribution in x-direction (N)");
legend("Case 1a","Case 1b", "Case 3");

subplot(1,2,2);
plot(span,M_z(1,:),"r", LineWidth = 1.5);
hold on
plot(span,M_z(2,:),"b--", LineWidth = 1.5);
plot(span,M_z(3,:),":", LineWidth = 1.5);
hold off
grid on
xlabel("Semi-span (m)");
ylabel("Bending moment distribution in z-direction (N)");
legend("Case 1a","Case 1b", "Case 3");


figure(6)
subplot(1,2,1);
plot(span,M_torque(1,:),"r", LineWidth = 1.5);
hold on
plot(span,M_torque(2,:),"b--", LineWidth = 1.5);
plot(span,M_torque(3,:),":", LineWidth = 1.5);
hold off
grid on
xlabel("Semi-span (m)");
ylabel(" (N)");
legend("Case 1a","Case 1b", "Case 3");

% Get the sectional load vs spanwise diagram for case 1&3 under different
% conditions so MTOW, EFW, MLW

%%%%%%% CHANGE THE VALUES OF MLW AND EFW
MTOW = 353385;         %mtow in kg
MLW =  320000        ; % maximum landing weoght
EFW =  300000      ; %EMPTY FUEL WEIGHT
n= 3.75;
mass=[MTOW, MLW, EFW];

for j= 1:3
    lift(j) = n* mass(j)*g;
    L0(j)= lift(j)/(32.5*pi*0.5 - 6.34);  %assume an elliptical lift distribution
    dl_new(j,:) = L0(j)*sqrt(1-(span/32.5).^2);
end

figure
plot(span, dl_new(1,:))
hold on
plot(span, dl_new(2,:))
plot(span, dl_new(3,:))
grid on
xlabel("Spanwise position (m)")
ylabel("Sectional load (N/m)")


figure
plot(span,dL(1,:)-dW_total(1,:),"r",LineWidth = 1.25);
hold on
plot(span,dL(3,:)-dW_total(3,:),LineWidth = 1.25,Color=[0 128 255]/255);
plot(span,dL(1,:),"--", LineWidth = 1.25,Color=[255 102 102]/255);
plot(span,dL(3,:),"--", LineWidth = 1.25,Color=[102 102 255]/255);
plot(span,-dW_total(1,:),":", LineWidth = 1.25,Color=[153 0 0]/255);
plot(span,-dW_total(3,:),":", LineWidth = 1.25,Color=[0 0 153]/255);
hold off
xlabel("Semi-span (m)");
ylabel("load distribution (N/m)");
grid on
legend("Case 1 - Resultant load","Case 3 - Resultant load","Case 1 - Lift ", "Case 3 - Lift ","Case 1 - Inertial loads","Case 3 - Inertial loads");


