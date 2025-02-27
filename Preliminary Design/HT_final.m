clear
clc
bluez=[0 0 1];
redz=[1 0 0];
greenz=[0 0.75 0];
lilaz=[0.75 0 0.75];
orangz=[1 0.5 0];
colorz={bluez, redz, greenz, orangz, lilaz};
n=1000; 
% Horizontal stabiliser planform dimensions:
S_ref = 482;
S=58; 
AR=5.8; 
lambda=0.4; 
bperp=sqrt(AR*S); % Span (perpendicular to fuselage)
sweep=35.155; 
b=bperp/cosd(sweep); % Length of wing
s=bperp/2; 
yspan = linspace(0,bperp/2,1000);% Getting spanwise points from root to tip
yroot=0; 
ytip=b/2;
save('yspan_htail', 'yspan')

xcg=36;
xcg_zf = 35.6;
xw=35.3;
% Getting chord distribution:
croot=(2*S)/(bperp*(1+lambda));
ctip=lambda*croot; 
c=((ctip-croot)/s).*yspan + croot;
xh=56; % Assuming ac is at 1/4 chord point

%% Getting total lift of horizontal stabilisers

% Getting lift from wings:
g=9.81; 
ULF=3.75; % Ultimate load factor (n=1.5*2.5=3.75)
MTOW=353000*g; 
MZFW = 162976*g; 
Lw=ULF*MTOW; 
Lw_zf = ULF*MZFW;

% Getting moment (lever) arms:
leverw=xcg-xw; 
leverh=xh-xcg; 

leverw_zf=xcg_zf-xw; 
leverh_zf=xh-xcg_zf; 

VD=214.06; % Dive speed (m/s, IAS)
VA = 171.25;
Cm0=-0.131; 
M0w=0.5*1.225*VD^2*S_ref*Cm0; % Zero lift pitching moment from wings
M0w_A=0.5*1.225*VA^2*S_ref*Cm0;

Lh=(Lw*leverw-M0w)/leverh; 
Lh_zf=(Lw_zf*leverw_zf-M0w)/leverh_zf; 

% SINCE SYMMETRICAL, ONLY NEED TO CONSIDER ONE WING:
Lh2=0.5*Lh; % Lift of one horizontal stabiliser in N
Lh2_zf=0.5*Lh_zf;

% Assuming elliptical distribution
LSecMax=(4*Lh2)/(pi*s); % in N
LSecMax_zf=(4*Lh2_zf)/(pi*s);

LSec=sqrt((LSecMax^2)*(1-(yspan.^2/s^2))); % Sectional lift at each spanwise point
LSec_zf=sqrt((LSecMax_zf^2)*(1-(yspan.^2/s^2)));

figure
box on
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1.5);

plot(yspan,LSec,'-b','LineWidth',1.5)
hold on
plot(yspan,LSec_zf,'-r','LineWidth',1.5)
xlabel('Spanwise Position [m]')
ylabel('Sectional Lift [N/m]')
legend('Lift (MTOW)','Lift (MZFW)')
grid on

hold off

lift=trapz(yspan,LSec)
lift2=trapz(yspan,LSec_zf)
%INERTIA

W=ULF*0.5*1230*g; % Weight of one horizontal stabiliser 

semispan = bperp/2;
S1 = 6.96;
S2 = S1*(lambda)^2;
WingBox_Volume = 1/3*(S1+S2+(S1*S2)^(1/2))*semispan;

Area_Span = (S1*(lambda^2-2*lambda+1))*(yspan./semispan).^2 + (S1*(-2+2*lambda)).*(yspan./semispan) + S1;
WingW_Span = ((W)/WingBox_Volume).*Area_Span;
wSec = WingW_Span;
wSecMax=(2*W)/(s*(ctip/croot + 1))
wSecMin=(ctip/croot)*wSecMax

figure
box on
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold','LineWidth', 1.5);

hold on
plot(yspan,LSec,'-b','LineWidth',1.5)
plot(yspan,LSec_zf,'-r','LineWidth',1.5)
plot(yspan,-WingW_Span,'-k',LineWidth=1.5)
xlabel('Spanwise Position [m]')
ylabel('Sectional Weight [N/m]')
grid on
legend('Lift (MTOW)','Lift (MZFW)','Weight')

hold off

dy=yspan(2);

for z=1:length(yspan) % z is spanwise coordinate but starting at the tip
    Shear(z) = dy*sum(LSec(length(LSec)-z+1:length(LSec))) - dy*sum(wSec(length(wSec)-z+1:length(wSec)));
    Shear_zf(z) = dy*sum(LSec_zf(length(LSec_zf)-z+1:length(LSec_zf))) - dy*sum(wSec(length(wSec)-z+1:length(wSec)));
end

% Flipping Shear to a coordinate system which starts from the root
Shear2=flip(Shear);
Shear2_zf=flip(Shear_zf);

figure
box on
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1.5);
hold on
plot(yspan,Shear2,'-b','LineWidth',1.5)
plot(yspan,Shear2_zf,'-r','LineWidth',1.5)
xlabel('Spanwise Position [m]')
ylabel('Shear Force [N]')
legend('MTOW Case','MZFW Case')
grid on
hold off



save('shear_force_htail.mat','Shear2','Shear2_zf','yspan')

%% Getting bending moment distribution of horizontal stabiliser 

for i=1:length(yspan) 
    BeamLift = dy*sum(LSec(length(LSec)-i+1:length(LSec))); 
    BeamLift_zf = dy*sum(LSec_zf(length(LSec_zf)-i+1:length(LSec_zf)));
    BeamWeight = dy*sum(wSec(length(wSec)-i+1:length(wSec))); 
    currentw = wSec(length(wSec)-i+1); 
    % Moment due to lift:
    MLift=(4/(3*pi))*BeamLift*(yspan(length(yspan))-yspan(length(yspan)-i+1));
    MLift_zf=(4/(3*pi))*BeamLift_zf*(yspan(length(yspan))-yspan(length(yspan)-i+1));
    % Moment due to weight:
    MWeight=((2*wSecMin+currentw)/(3*(wSecMin +currentw)))*BeamWeight*(yspan(length(yspan))-yspan(length(yspan)-i+1));
    % Bending Moment:
    BM(i)=MWeight - MLift;
    BM_zf(i)=MWeight - MLift_zf;
end

% Flipping BM to a coordinate system which starts from the root:
BM2=flip(BM);


BM2_zf=flip(BM_zf);



figure
box on
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1.5);
hold on
grid on
plot(yspan,BM2,'-b','LineWidth',1.5)
plot(yspan,BM2_zf,'-r','LineWidth',1.5)
xlabel('Spanwise Position [m]')
ylabel('Bending Moment [Nm]')
legend('MTOW Case','MZFW Case')
hold off

save("bm_htail", "BM2_zf")

%% Getting torque distribution of horizontal stabiliser

AC=0.25; 
FA=0.45; % Flexural analysis
CG=0.391809984147464; 

% Finding torque distribution:
for i=1:length(yspan)
    SectionLiftMoment(i)=LSec(length(LSec)-i+1)*((FA-AC)*c(length(c)-i+1)); 
    SectionLiftMoment_zf(i)=LSec_zf(length(LSec)-i+1)*((FA-AC)*c(length(c)-i+1));
    SectionWeightMoment(i)=wSec(length(wSec)-i+1)*((CG-FA)*c(length(c)-i+1));
    TorqueSec(i)=SectionLiftMoment(i)+SectionWeightMoment(i); 
    Torque(i)= dy*sum(TorqueSec(length(TorqueSec)-i+1:length(TorqueSec))); 
    TorqueSec_zf(i)=SectionLiftMoment_zf(i)+SectionWeightMoment(i); 
    Torque_zf(i)= dy*sum(TorqueSec_zf(length(TorqueSec_zf)-i+1:length(TorqueSec_zf))); 
end

Torque2=flip(Torque);
Torque2_zf=flip(Torque_zf);


figure
box on
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1.5);
hold on
grid on
plot(yspan,Torque2,'-b','LineWidth',1.5)
plot(yspan,Torque2_zf,'-r','LineWidth',1.5)
xlabel('Spanwise Position [m]')
ylabel('Torque [Nm]')
ylim([0,200000])
xlim([0,10])
legend('MTOW Case','MZFW Case')
hold off

torqueData = Torque2;                     
save("torque_htail.mat", "torqueData"); 
maxvalue = [max(abs([Shear,Shear2_zf])),max(abs([BM2,BM2_zf])),max(abs([Torque,Torque_zf]))]

%wing box geometry & load estimation input

yspan = load("yspan_htail.mat"); %array of y coords spanwise from load code
yspan = yspan.yspan;
bm = load("bm_htail.mat"); %array of bending moment from root to tip in Nm from load code
bm = bm.BM2_zf;

%material properties (whatever the material is here!) 
E_r = 84; %Young's modulus (GPa) of ribs
E = 74.8e9; %Young's modulus (Pa) of skin
sigma_y_r = 450; %yield stress (MPa) of ribs
sigma_y = 756; %yield stress (MPa) of skin-stringer

n_ribs = 20;
b_h =14.1/2;
L = b_h/n_ribs %Rib spacing in m

%Wing parameters
%wing parameter

%wing chord distribution starting at the wing-fuselage intersection
syms y
c(y) = 7.5 + (3.0 - 7.5)*(y/7.05); %for our tailplane HT

% defining span stations
s_station = linspace(0, b_h, n_ribs);

% defining new bm matching the span stations
bm_new = abs(spline(yspan, bm, s_station));

c_box = double(vpa(0.5 * c(s_station))); %wingbox width (m)
b_box = double(vpa(0.0942 * c(s_station))); %wingbox average height (m)
N = bm_new./(c_box.*b_box); %compressive stress in N/m


%1. Initial buckling
t_init = ((N .* b_box.^2)./(3.62 * E)).^(1/3) * 1000; 
cr_stress = N./(1000*t_init); 
area_init = c_box .* t_init .* 1000; 


%2. Design panel size and stringer size
n_st = 22;
W = [];
for j = 1:length(n_st)
%initialising the input values
n = n_st(j); %number of panels
ts = 1.5; %stringer thickness in mm
h = 32; %stringer web height in mm
f2w = 0.3; %stringer flange to web ratio (this is typically 0.3 in the lectures, but higher values exist such as 0.4 and 0.5 in the ESDU database.)

%calculations from the inputs
d = f2w * h; %flange width of stringer
As = (h + (2*d)) * ts; %cross-sectional area of the stringer
b = c_box/n * 1000; %panel width in mm

%3. Use initial buckling of panel to find skin thickness
t = 1000 .* ((N .* (b/1000).^2) ./ (3.62 * E)) .^ (1/3); %skin thickness with stringers added in mm
cr_s0 = N./(1000*t); %critical buckling stress in N/mm^2 (sigma0)
t_e = t + As./b; %panel effective thickness in mm (smeared thickness?)
As_e = b .* n .* t_e; %panel effective area in mm^2

%R_t = ts./t; %ratio of ts/t
R_t = ts./t(1)
%R_b = As./b./t; %ratio of As/bt
R_b = As./(b(1).*t(1))


%--- Calculate the ratios ---
% Use the first element of t and b for the current design point
R_t = ts ./ t(1);   % ratio of stringer thickness to skin thickness
R_b = As ./ (b(1) * t(1));  % ratio of stringer area to (panel width * skin thickness)

% Clamp R_t to the valid range for the Catchpole diagram
R_t = max(min(R_t, 1.4), 0.4);
% Optionally, ensure R_b is within the range [0, 1.2]:
R_b = max(min(R_b, 1.2), 0);

tableCF = readtable("catchpole_data.csv");

catchpoleBlock = tableCF(1:8, 1:10);
farrarBlock    = tableCF(9:16, 1:10);

catchpoledata_matrix = table2array(catchpoleBlock);
farrardata_matrix    = table2array(farrarBlock);

stressFactors = catchpoledata_matrix(2:end, 2:end);  % for Catchpole
f_testing     = farrardata_matrix(2:end, 2:end);        % for Farrar

X_grid = [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2];             % 7 points (for R_b)
Y_grid = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4];   % 9 points (for R_t)

%% Interpolate the diagram values
% Transpose the matrices so that they are 9×7, matching [length(Y_grid) x length(X_grid)]
R_sigma = interp2(X_grid, Y_grid, stressFactors', R_b, R_t);
F       = interp2(X_grid, Y_grid, f_testing',    R_b, R_t);

% Compute the adjusted critical buckling stress:
cr_stress1 = R_sigma .* cr_s0;


%5. Getting FARRAR efficiency for each span station

mean_stress = F.*sqrt((N .* E*10^9)./ L) * 10^-6;
L = N(1) * E / (cr_stress1(1) * 10^6/F)^2 % length in m

end

%7. Plotting for skin thickness and stress distribution
t_step = ceil(t*2)/2;

figure
hold on
plot(s_station, t, '--r', 'LineWidth', 1.5)
stairs(s_station ,t_step, 'b', 'LineWidth', 1.5)
grid minor
xlabel("Span length [m]")
ylabel("Skin thickness [mm]")
xlabel(get(get(gca,'XLabel'),'String'),'Interpreter','latex');
ylabel(get(get(gca,'YLabel'),'String'),'Interpreter','latex');
lg = legend("Theoretical skin thickness [mm]", "Manufactured skin thickness [mm]");
% Set global properties for all axes
set(findobj(gcf, 'type', 'axes'),'FontSize', 13, 'FontWeight', 'Bold', 'LineWidth', 1);
set(findobj(gcf, 'type', 'line'), 'LineWidth', 1.5);
set(lg, 'Interpreter','latex');
xlim([0 b_h-0.01])
ylim([0 4])
hold off


%8. Interaction between shear and compression
yspan = linspace(0,b_h, 1000);

% 1) Load your torque .mat file
torqueData = load("torque_htail.mat");  % or whatever it's called
yspanTorque = yspan;
Tval        = torqueData.torqueData;         % torque array [N·m]

% 2) Interpolate torque to your stations
T_s = spline(yspanTorque, Tval, s_station);

% 3) Compute local q0 at each station
q0 = T_s ./ (2 .* c_box .* b_box .* 1000);   % [N/mm]

% 4) Shear stress from torque
s_q0 = q0 ./ t;  % [N/mm^2]

tresca_yield = 378

R_s = s_q0 ./ tresca_yield
R_c = cr_s0 ./ cr_stress1

check = R_s .^2 + R_c

figure
set(findobj(gcf,'type','axes'),'FontSize',12,'FontWeight','Bold', 'LineWidth', 1.5);
plot(s_station, check, 'b', 'LineWidth', 1.5 )
grid minor
xlabel("Span length (m)", Interpreter="latex")
ylabel("$R_s^2 + R_c$ (mm)", Interpreter="latex")
xlim([0 5.2])
% Set global properties for all axes
set(findobj(gcf, 'type', 'axes'),'FontSize', 13, 'FontWeight', 'Bold', 'LineWidth', 1);
set(findobj(gcf, 'type', 'line'), 'LineWidth', 1.5);
set(lg, 'Interpreter','latex');
xlabel(get(get(gca,'XLabel'),'String'),'Interpreter','latex');
ylabel(get(get(gca,'YLabel'),'String'),'Interpreter','latex');
box on
% Skin-stringer panel weight
rho = 2850; % density of AL2024 in kg/m^3
skin_stringer_weight = sum(As_e .* L * rho / 1e6); % weight in kg

% Load Shear Force and Torque Data for the horizontal tail
shearData_ht = load("shear_force_htail.mat"); % Shear force along span
torqueData_ht = load("torque_htail.mat");     % Torque along span

% Display available field names (optional check)
disp(fieldnames(shearData_ht));
disp(fieldnames(torqueData_ht));

% Extract shear and torque data
shear_span_ht = shearData_ht.yspan;
shear_force_ht = shearData_ht.Shear2;


% Interpolate shear force and torque to match horizontal tail span stations
V_new_ht = spline(shear_span_ht, shear_force_ht, s_station); 
T_new_ht = T_s

% Define Material Properties for Spar Webs
E_Aluminium = 72.5;   % GPa (Al2024)
Ks_spar     = 8.5;    % Buckling coefficient for spars

% Preallocate Arrays for Web Thickness at Front and Rear Spars
t_FS_ht = zeros(size(s_station)); % Front Spar Theoretical
t_RS_ht = zeros(size(s_station)); % Rear Spar Theoretical

% Compute Spar Web Thickness using Buckling Formula
for i = 1:length(s_station)
    % Shear flows [N/mm] for front and rear spars
    q0_ht = T_new_ht(i) / (2 * c_box(i) * b_box(i) * 1000); 
    q2_ht = V_new_ht(i) / (2 * b_box(i) * 1000);

    q_FS_ht = q2_ht + q0_ht;  % Front spar
    q_RS_ht = q2_ht - q0_ht;  % Rear spar

    % Spar Web Thickness Calculation from Buckling Formula (in mm)
    t_FS_ht(i) = (q_FS_ht * 1000 * b_box(i) / (Ks_spar * E_Aluminium * 1e9))^(1/3) * 1000;
    t_RS_ht(i) = (q_RS_ht * 1000 * b_box(i) / (Ks_spar * E_Aluminium * 1e9))^(1/3) * 1000;
end

% Define Manufactured Thickness (Rounded to Nearest 0.5mm)
t_FS_manu = ceil(t_FS_ht * 2) / 2;  % Manufactured Front Spar
t_RS_manu = ceil(t_RS_ht * 2) / 2;  % Manufactured Rear Spar
spar_weight = sum((t_FS_manu + t_RS_manu) .* b_box .* L * rho / 1e6); % weight in kg
%% Plot Spar Web Thickness Distribution for Horizontal Tailplane
figure
hold on

% Theoretical Thickness (Smooth Lines)
plot(s_station, t_FS_ht, '--b', 'LineWidth', 1.5, 'DisplayName','Theoretical Front Spar')
plot(s_station, t_RS_ht, '--r', 'LineWidth', 1.5, 'DisplayName','Theoretical Rear Spar')

% Manufactured Thickness (Step Plot)
stairs(s_station, t_FS_manu, '-b', 'LineWidth', 2, 'DisplayName','Manufactured Front Spar')
stairs(s_station, t_RS_manu, '-r', 'LineWidth', 2, 'DisplayName','Manufactured Rear Spar')

hold off
grid minor
xlabel("Spanwise Location (m)", 'Interpreter','latex', 'FontSize', 14)
ylabel("Spar Web Thickness (mm)", 'Interpreter','latex', 'FontSize', 14)
legend('Location','best','Interpreter','latex', 'FontSize', 12)
title("Spar Web Thickness Distribution - Horizontal Tailplane",'Interpreter','latex', 'FontSize', 14)

disp("=== Spar Web Thickness Calculation for H-Tail Completed ===");

%Rib Design
% iterating through number of ribs for the weight optimisation
%n_ribs = 15:1:30;
% rib_array = 0.05:0.05:0.7;
rib_array = 0.44;


for i = 1 : length(rib_array)
L = rib_array(i);
% no. of ribs
n_ribs = b_h/L;

% defining span stations
s_station = linspace(0,b_h, n_ribs);

%wing chord distribution starting at the wing-fuselage intersection
syms y
c(y) = 7.5 + (3.0 - 7.5)*(y/7.05);
chord = double(vpa(c(s_station)));

% defining new bm matching the span stations
bm_new = spline(yspan, bm, s_station);

%solve this problem of bm_new at the end going to the negative which causes
%problems of having complex number later on

c_box = double(vpa(0.5 .* chord)); %wingbox width
b_box = double(vpa(0.0942 .* chord)); %wingbox average height (m)
N = bm_new./(c_box.*b_box); %compressive stress in N/m

%1. Rib thickness
% calculating effective skin thickness
n = 16; %number of panels
ts = 1.5; %stringer thickness in mm
h = 32; %stringer web height in mm
f2w = 0.3; %stringer flange to web ratio (this is typically 0.3 in the lectures, but higher values exist such as 0.4 and 0.5 in the ESDU database.)

%calculations from the inputs
FARRAR = 0.75;
d = f2w * h; %flange width of stringer
As = (h + (2*d)) * ts; %cross-sectional area of the stringer
b = c_box/n * 1000; %panel width in mm
t = 1000 .* ((N .* (b/1000).^2) ./ (3.62 * E)) .^ (1/3); %skin thickness with stringers added in mm
t_e = t + As./b; %panel effective thickness in mm (smeared thickness?)
As_e = b .* n .* t_e; %panel effective area in mm^2

%% Calculating the second moment of area of skin panel

Is = ((c_box.*(t_e./1000).^3)/12) + (c_box.*(t_e./1000).*(b_box./2).^2);

% Calculating crush load

crush_load = ((bm_new.^2) .* L .* b_box .* (t_e./1000) .* c_box) ./ (2 .* E_r*10^9 .* (Is.^2)); 

% rib thickness from buckling stress
t_rib_b = ((crush_load .* b_box.^2) ./ (3.62 * E_r * 10^9 .* c_box)) .^(1/3) .* 1000 % in mm

% rib thickness from yield stress
t_rib_y = (crush_load ./ ((sigma_y_r*10^6) .* c_box)) * 1000; %(mm)

% crush stress in 
sigma_crit(i) = max(FARRAR .* sqrt((N .* E) ./ L));

%t_smeared_total = t_e + (((b_box./1000) .* t_rib_b) ./ (L / 1000));

rho = 2765; %density of skin material in kg/m^3
rho_r = 2765; %density of rib material in kg/m^3
rib_weight = sum(t_rib_b .* b_box .* L * rho_r / 1e6); % weight in kg

end



figure
hold on
plot(s_station, t_rib_b,'xB',"LineWidth",1.5, 'MarkerSize',6)
plot(s_station, ones(1, length(s_station)), 'xr')
ylim([0 1.3])
xlim([0 b_h])
grid minor
ylabel("Rib Thickness (mm)",'Interpreter','latex')
xlabel("Spanwise location (m)",'Interpreter','latex')
lg = legend(' Theoretical thickness (mm)', 'Manufactured thickness (mm)')
% Set global properties for all axes
set(findobj(gcf, 'type', 'axes'),'FontSize', 13, 'FontWeight', 'Bold', 'LineWidth', 1);
set(findobj(gcf, 'type', 'line'), 'LineWidth', 1.5);
xlabel(get(get(gca,'XLabel'),'String'),'Interpreter','latex');
ylabel(get(get(gca,'YLabel'),'String'),'Interpreter','latex');
set(lg, 'Interpreter','latex');

box on
hold off

%D-section design
N = 20; % number of pseudo-ribs;
t_guess = 1 ;% mm
t1 = t_guess.*ones(1,N); % guess d cell thickness array , mm

%Shear calculation
semi_span =  b_h; % wing semi span outside fuselage, m
bref = 0:(semi_span/N):semi_span - (semi_span/N);
E =  79.8*10^3; %Young's modulus, MPa

box_h = b_box(1) ; % wing box height, m
tresca_stress = 320 % yield shear stress for Al 2090, MPa

a = semi_span/(N+1) % d cell span length, m

%wing chord distribution starting at the wing-fuselage intersection
syms y
c(y) = 7.5 + (3.0 - 7.5)*(y/7.05); %for our tailplane HT

b = double(0.15 .* c(bref)); % width of d cell, m


R = double(0.5.*((b.^2)./(box_h/2) + (box_h/2)^2./b)); % average radius of curvatures array, m

a/b(1)
b(1)/a

coeff_1 = b./sqrt(R(1).*(t1./1000))  % for a > b
coeff_2 = a./sqrt(R(1).*(t1./1000)) % for a < b


%Getting Ks from plot
x = [0, 2, 4, 6 ,8 ,10, 12, 14];
y = [9, 12, 17, 23, 29, 37, 45, 54];

Ks = abs(spline(x,y,coeff_1))


tau_crit = Ks(1).*E.*(t1(1)./(b(1).*1000)).^2 % critical buckling stress array, MPa, we want tau_crit = tresca_stress

t_req = sqrt(tresca_stress./(Ks.*E)).*b*1000

%Plotting
figure
set(findobj(gcf,'type','axes'),'FontSize',20,'FontWeight','Bold', 'LineWidth', 1.5);
plot(bref, t_req, 'xk', LineWidth=1.5)
grid minor
xlabel('Span length (m)', Interpreter="latex", FontSize=14)
ylabel("Pseudo-ribs thickness (mm)", Interpreter="latex", FontSize=12)
xlim([0 7.05])
box on


%WING PLOT
semispan_vt = 14.1/2;
n = 500;
rib_no = 14;
stringer_no = 17;

lambda = 0.4;
sweep_v_quarter = 31.6;
wing_root_vt = 6.4;
wing_tip_vt = 0.4*6.8;
sweep_v_LE = atand(tand(sweep_v_quarter) - (4 / 1.6) * ((0 - 25) / 100) * ((1 - lambda) / (1 + lambda)));
fspar_pos = 0.2;
rspar_pos = 0.7;

rib_pos = linspace(0,semispan_vt,rib_no);

yspan_vt = linspace(0,semispan_vt,n);

LE_x = yspan_vt./tand(90-sweep_v_LE);

root_x = linspace(0,wing_root_vt,n);
root_y = zeros([1,n]);

tip_x = (semispan_vt/tand(90-sweep_v_LE))+linspace(0,wing_tip_vt,n);
tip_y = semispan_vt.*ones([1,n]);

TE_x = linspace(wing_root_vt,max(tip_x),n);

F_spar_x = linspace(fspar_pos*wing_root_vt,(semispan_vt/tand(90-sweep_v_LE))+fspar_pos*wing_tip_vt,n);
R_spar_x = linspace(rspar_pos*wing_root_vt,(semispan_vt/tand(90-sweep_v_LE))+rspar_pos*wing_tip_vt,n);

stringer_pos = linspace(fspar_pos*wing_root_vt,rspar_pos*wing_root_vt,stringer_no);
stringer_pos_tip = linspace(max(LE_x)+fspar_pos*wing_tip_vt,max(LE_x)+rspar_pos*wing_tip_vt,stringer_no);

figure 

%WING PLOTS:
plot(root_x,root_y,'k',LineWidth = 1.0)
hold on
grid on
plot(tip_x,tip_y,'k',LineWidth = 1.0)
plot(LE_x,yspan_vt,'k',LineWidth = 1.0)
plot(TE_x,yspan_vt,'k',LineWidth = 1.0)

%SPARS
plot(F_spar_x,yspan_vt,'Color','#FF0000',LineWidth = 2.5)
plot(R_spar_x,yspan_vt,'Color','#FF0000',LineWidth = 2.5)

%RIBS
for i = 2:rib_no-1
    line([rib_pos(i)/(tand(90-sweep_v_LE)),rib_pos(i)/(tand(90-sweep_v_LE))+rspar_pos*(wing_root_vt+rib_pos(i)*(wing_root_vt*(lambda-1))/(semispan_vt))],[rib_pos(i),rib_pos(i)],'Color','#0000FF','LineWidth',1.8)
end

%STRINGERS
for i = 2:stringer_no-1
    line([stringer_pos(i),stringer_pos_tip(i)],[0,semispan_vt],'Color','#00FF00','LineWidth',0.1)
end
h = zeros(3, 1);
h(1) = plot(NaN,NaN,Color = '#FF0000');
h(2) = plot(NaN,NaN,Color = '#0000FF');
h(3) = plot(NaN,NaN,Color = '#00FF00');
legend(h, 'H-Tail Spars','H-Tail Ribs','H-tail Stringers','Location','Northwest');

ylim([0,10])
set(findobj(gcf, 'type', 'axes'),'FontSize', 13, 'FontWeight', 'Bold', 'LineWidth', 0.05);
set(findobj(gcf, 'type', 'line'), 'LineWidth', 1.5);
xlabel(get(get(gca,'XLabel'),'String'),'Interpreter','latex');
ylabel(get(get(gca,'YLabel'),'String'),'Interpreter','latex');
lg = legend;
set(lg, 'Interpreter','latex');
ylim([0,8])
axis equal
grid on
xlabel('X Position (m)')
ylabel('Y Position (m)')
hold off
