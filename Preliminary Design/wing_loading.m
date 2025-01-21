% This code calculates the loads acting on the wing

%% Housekeeping

clear
clc
close all


%% Defining some key values (wing only)

load("parameters.mat");

% Case
MTOW = 45530; % kg
MZFW = 0; % kg
Weight = MTOW; % MTOW or MZFW
lf = 3.75; % load factor

%masses
wdrymass = 3468.5; %dry mass of the wing
wfuelmass = 4522.1; %wing fuel mass
fuelden = 775; %fuel density
wfueltankl = 10.67; %wing fuel tank length per wing
emass = 1620; %mass of the engine
maingearmass = 456; %mass of main gear

%aerodynamics (assuming an elliptic lift distribution)
Lrequired = Weight*lf*9.81;
L0 = (2*Lrequired)/(pi*12.42); %lift at centre of ellipse (calculated in excel)

%aerofoil properties
cog = 0.385; %location of cog of the aerofoil (in terms of chord lengths)
larm = 0.135; %location of cog from AC (in terms of chord lengths)
frontspar=0.1; %location of front spar from LE (in terms of chord lengths)
rearspar=0.65; %location of rear spar from TE (in terms of chord lengths)
flexax=(frontspar+rearspar)/2; %location of flexural axis (shear centre) (in terms of chord lengths)

%wing properties
syms x
c(x) = 4.9469-0.2610*(x+1.6);

%other stuff
spanwofuse = param.wing.b/2 - 3.2/2; %span of wing (outside fuselage)
g = 9.81; %acceleration due to gravity

%% Dicretise the semi-span of the wing

spanstation = [0:0.05:spanwofuse];

%% Wing inertia load

%calculating loads on the wing (with a load factor of 2.5)
wload = (wdrymass./spanwofuse)*lf*g; %load due to mass of the wing
fload = wfuelmass*g*lf/wfueltankl; %load due to fuel in the wing
eload = emass*g*lf; %load due to engine
mgload = maingearmass*g*lf; %load due to main gear

%put each of the loads in the right dicretisation in the wing
for i = 1:length(spanstation)
    w(i) = wload;

    if (spanstation(i)>1.8478) && (spanstation(i)<11.9245)
        fuelload(i) = fload;
        if (spanstation(i) == 2.6)
            engineload(i) = eload;
        end
    elseif (spanstation(i) == 0.4)
        maingload(i) = mgload;
    else
        fuelload(i) = 0;
        engineload(i) = 0;
        maingload(i) = 0;
    end
end

%calculating the total load (at each section) due to the inertia on the wing
for i = 1:length(spanstation)-1
    totload(i) = ((w(i+1) + w(i))*(spanstation(i+1)-spanstation(i))/2) + ((fuelload(i+1) + fuelload(i))*(spanstation(i+1)-spanstation(i))/2) + engineload(i) + maingload(i);
end
totload(length(spanstation)) = 0;

%calculating shear due to loads on each section of the wing
for i = 1:length(spanstation)
    temp = totload(i:length(spanstation));
    shear(i) = sum(temp);
end

%calculating sectional moment due to shear
for i = 1:length(spanstation)-1
    dM(i) = (shear(i+1)+shear(i))*(spanstation(i+1)-spanstation(i))/2;
end
dM(length(spanstation)) = 0;

%calculating total moment due to sectional moment
for i = 1:length(spanstation)
    temp = dM(i:length(spanstation));
    bendmom(i) = sum(temp);
end

inertia = [spanstation',w',fuelload',engineload',maingload',totload',shear',dM',bendmom'];



%plot inertia stuff

% % 
% %plot shear force distribution
% figure
% plot(spanstation,shear)
% legend('Shear force','Interpreter','latex')
% xlabel('Spanwise distance from root, m','Interpreter','latex')
% ylabel('Shear force, N','Interpreter','latex')
% grid
% box on
% ax = gca;
% ax.FontSize = 15;
% set(gcf,'units','inches','position',[1,1,8,6])
% % 
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% 
% %plot bending moment distribution
% figure
% plot(spanstation,bendmom)
% legend('Bending moment','Interpreter','latex')
% xlabel('Spanwise distance from root, m','Interpreter','latex')
% ylabel('Bending moment, Nm','Interpreter','latex')
% grid
% box on
% ax = gca;
% ax.FontSize = 15;
% set(gcf,'units','inches','position',[1,1,8,6])
% % 
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];

%% Wing Aerodynamic Loads

%% Symmetric load case

%calculating sectional lift (N/m)
for i = 1:length(spanstation)
    dL(i) = L0*sqrt(1-(spanstation(i)/spanwofuse)^2);
end

%calculating lift
for i = 1:length(spanstation)-1
    L(i) = (dL(i+1)+dL(i))*(spanstation(i+1)-spanstation(i))/2;
end
L(length(spanstation)) = 0;

%calculating shear force (due to aerodynamic loads)
for i = 1:length(spanstation)
    temp = L(i:length(spanstation));
    shearae(i) = sum(temp);
end

%calculating sectional moment from shear (due to aerodynamic loads)
for i = 1:length(spanstation)-1
    dMae(i) = (shearae(i+1)+shearae(i))*(spanstation(i+1)-spanstation(i))/2;
end
dMae(length(spanstation)) = 0;

%calculating bending moment (due to aerodynamic loads)
for i = 1:length(spanstation)
    temp = dMae(i:length(spanstation));
    bendmomae(i) = sum(temp);
end

symaeload = [spanstation',dL',L',shearae',dMae',bendmomae'];

%plotting stuff

% %lift distribution
% figure
% plot(spanstation,dL)
% 
% %shear force
% figure
% plot(spanstation,shearae)
% 
% %bending moment
% figure
% plot(spanstation,bendmomae)


%% Aero + Inertia Loading on the wing

wingload = [spanstation', shear'+shearae', bendmom'+bendmomae'];

%% Wing torque

%find the chord distribution (m)
for i = 1:length(spanstation)
    cdis(i) = double(vpa(c(spanstation(i))));
end

%moment due to sectional lift (from cog) (Nm/m)
aerola = dL*larm.*(cdis - (tand(20.15)).*(spanwofuse - spanstation));

%conditions for the first segment (heaviest)
[T,a,p,rho] = atmosisa(convlength(25000,'ft','m'));

%torque due to zero lift moment coefficient
M0 = -0.079*(0.5*rho*143^2)*cdis*param.wing.S; %idk if the C_m0 is correct (CHECK WITH XFOIL)

%torque due to wing weight
wtorque = (cog-flexax)*cdis.*wload;

%torque due to fuel weight
cogc = cog*cdis; %location of cog of the wing (for each spansection)
fcog(x) = vpa(3.5865145567502431100948001585493 - 0.36047171743782140604372976611103*(x+1.6)); %COG of the fuel (from the LE)

fuelt = (double(fcog(cdis)) - cogc).*fuelload; %moment due to fuel load (about COG of wing)

%torque due to engine
ecog = 0; %location of cog of engine on the wing
engt = (ecog-cogc).*engineload;

%torque due to undercarriage
ucog = 0.55; %location of cog of uc on the wing (in terms of c)
ut = (ucog-cogc).*maingload;


% finding the difference in moments between two span stations and hence the
% sectional torque
for i = 1:length(spanstation)-1
    aerola1(i) = (aerola(i)+aerola(i+1))*(spanstation(i+1)-spanstation(i))/2;
    M01(i) = (M0(i)+M0(i+1))*(spanstation(i+1)-spanstation(i))/2;
    wtorque1 = (wtorque(i)+wtorque(i+1))*(spanstation(i+1)-spanstation(i))/2;
    fuelt1(i) = (fuelt(i)+fuelt(i+1))*(spanstation(i+1)-spanstation(i))/2;
    engt1(i) = (engt(i)+engt(i+1))*(spanstation(i+1)-spanstation(i))/2;
    ut1(i) = (ut(i)+ut(i+1))*(spanstation(i+1)-spanstation(i))/2;
end
aerola1(length(spanstation)) = 0;
M01(length(spanstation)) = 0;
wtorque1(length(spanstation)) = 0;
fuelt1(length(spanstation)) = 0;
engt1(length(spanstation)) = 0;
ut1(length(spanstation)) = 0;

%torque distribution
for i = 1:length(spanstation)
    combT(i) = aerola1(i) + M01(i) + wtorque1(i) + fuelt1(i) + engt1(i) + ut1(i);
end

for i = 1:length(spanstation)
    temp = combT(i:length(spanstation));
    torque(i) = sum(temp);
end
torque(length(spanstation)) = 0;

%plotting stuff

% figure
% plot(spanstation,combT)
% grid minor
% title("Torque distribution")
% 
figure
plot(spanstation,torque)
grid minor
title("Torque diagram")

wingtorque = [spanstation', combT', torque'];


%% Nose off case
% Load factor is 3, no aerodynamic loads, point load acting at main gear
% location of the wing (3*MTOW)

% Case
MTOW = 45530; % kg
MZFW = 0; % kg
Weight = MTOW; % MTOW or MZFW
lf = 3; % load factor

fload = wfuelmass*g*lf/wfueltankl; %load due to fuel in the wing
eload = emass*g*lf; %load due to engine
mgload = maingearmass*g*lf; %load due to main gear

%% Wing inertia load

%calculating loads on the wing (with a load factor of 3)
wload = (wdrymass./spanwofuse)*lf*g; %load due to mass of the wing
fload = wfuelmass*g*lf/wfueltankl; %load due to fuel in the wing
eload = emass*g*lf; %load due to engine
mgload = maingearmass*g*lf + 3*MTOW*g; %load due to main gear

%put each of the loads in the right dicretisation in the wing
for i = 1:length(spanstation)
    w1(i) = wload;
    if (spanstation(i)>1.8478) && (spanstation(i)<11.9245)
        fuelload(i) = fload;
        if (spanstation(i) == 2.6)
            engineload(i) = eload;
        end
    elseif (spanstation(i) == 0.4)
        maingload(i) = mgload;
    else
        fuelload(i) = 0;
        engineload(i) = 0;
        maingload(i) = 0;
    end
end

%calculating the total load (at each section) due to the inertia on the wing
for i = 1:length(spanstation)-1
    totload(i) = ((w1(i+1) + w1(i))*(spanstation(i+1)-spanstation(i))/2) + ((fuelload(i+1) + fuelload(i))*(spanstation(i+1)-spanstation(i))/2) + engineload(i) + maingload(i);
end
totload(length(spanstation)) = 0;


%calculating shear due to loads on each section of the wing
for i = 1:length(spanstation)
    temp = totload(i:length(spanstation));
    shear(i) = sum(temp);
end

%calculating sectional moment due to shear
for i = 1:length(spanstation)-1
    dM(i) = (shear(i+1)+shear(i))*(spanstation(i+1)-spanstation(i))/2;
end
dM(length(spanstation)) = 0;

%calculating total moment due to sectional moment
for i = 1:length(spanstation)
    temp = dM(i:length(spanstation));
    bendmom(i) = sum(temp);
end

inertia1 = [spanstation',w1',fuelload',engineload',maingload',totload',shear',dM',bendmom'];

wingload1 = [spanstation', shear', bendmom'];

%% Wing torque

%find the chord distribution (m)
for i = 1:length(spanstation)
    cdis(i) = double(vpa(c(spanstation(i))));
end


%conditions for the first segment (heaviest)
[T,a,p,rho] = atmosisa(convlength(25000,'ft','m'));

%torque due to zero lift moment coefficient
M0 = -0.079*(0.5*rho*143^2)*cdis*param.wing.S; %idk if the C_m0 is correct (CHECK WITH XFOIL)

%torque due to wing weight
wtorque = (cog-flexax)*cdis.*wload;

%torque due to fuel weight
cogc = cog*cdis; %location of cog of the wing (for each spansection)
fcog(x) = vpa(3.5865145567502431100948001585493 - 0.36047171743782140604372976611103*(x+1.6)); %COG of the fuel (from the LE)

fuelt = (double(fcog(cdis)) - cogc).*fuelload; %moment due to fuel load (about COG of wing)

%torque due to engine
ecog = 0; %location of cog of engine on the wing
engt = (ecog-cogc).*engineload;

%torque due to undercarriage
ucog = 0.55; %location of cog of uc on the wing (in terms of c)
ut = (ucog-cogc).*maingload;


% finding the difference in moments between two span stations and hence the
% sectional torque
for i = 1:length(spanstation)-1
    M01(i) = (M0(i)+M0(i+1))*(spanstation(i+1)-spanstation(i))/2;
    wtorque1 = (wtorque(i)+wtorque(i+1))*(spanstation(i+1)-spanstation(i))/2;
    fuelt1(i) = (fuelt(i)+fuelt(i+1))*(spanstation(i+1)-spanstation(i))/2;
    engt1(i) = (engt(i)+engt(i+1))*(spanstation(i+1)-spanstation(i))/2;
    ut1(i) = (ut(i)+ut(i+1))*(spanstation(i+1)-spanstation(i))/2;
end
M01(length(spanstation)) = 0;
wtorque1(length(spanstation)) = 0;
fuelt1(length(spanstation)) = 0;
engt1(length(spanstation)) = 0;
ut1(length(spanstation)) = 0;

%torque distribution
for i = 1:length(spanstation)
    combT1(i) = M01(i) + wtorque1(i) + fuelt1(i) + engt1(i) + ut1(i);
end

for i = 1:length(spanstation)
    temp = combT1(i:length(spanstation));
    torque1(i) = sum(temp);
end
torque1(length(spanstation)) = 0;

wingtorque1 = [spanstation', combT1', torque1'];

%% Upload to the main folder to be used in other codes

save('wingload1.mat', 'wingload') %for symmetric load case
save('wingtorque1.mat', 'wingtorque') %for symmetric load case

save('wingload2.mat', 'wingload1') %for nose off load case
save('wingtorque2.mat', 'wingtorque1') %for nose off load case

% 
% figure 
% hold on
% plot(spanstation, wingload(:,3))
% plot(spanstation, wingload1(:,3))
% legend("wingload", "wingload1")


plot(spanstation,wingload(:,3))
legend('Bending moment','Interpreter','latex')
xlabel('Spanwise distance from root, m','Interpreter','latex')
ylabel('Bending moment, Nm','Interpreter','latex')
grid
box on
ax = gca;
ax.FontSize = 15;
set(gcf,'units','inches','position',[1,1,8,6])
% 
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];



figure
hold on
plot(spanstation,L)
plot(spanstation,-w-fuelload)
plot(spanstation,-w-fuelload+L)
hold off
legend("Aerodynamic Loads","Inertial Loads (excluding landing gear and engine)","Total Load",'Interpreter','latex')
xlabel('Spanwise distance from root, m','Interpreter','latex')
ylabel('Load, N','Interpreter','latex')
ylim([-2.5*10^4 1.2*10^4])

grid
box on
ax = gca;
ax.FontSize = 15;
set(gcf,'units','inches','position',[1,1,8,6])
% 
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
