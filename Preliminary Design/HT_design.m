% housekeeping
%catchpole.m and farraval.m will change based on the composite material
%chosen we can use a graph reader for catchpole to get x and y values and
%the other is using farrar efficiency 
clear,clc,close all

%Skin-stringer
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

n_ribs = 12;
b_h =11.32/2;
L = b_h/n_ribs %Rib spacing in m

Wing parameters
%wing parameter

%wing chord distribution starting at the wing-fuselage intersection
syms y
c(y) = 7.5 + (3.0 - 7.5)*(y/7.05); %for our tailplane HT

% defining span stations
s_station = linspace(0, b_h, n_ribs);

% defining new bm matching the span stations
bm_new = spline(yspan, bm, s_station);
bm_new(bm_new <0)= 0;

c_box = double(vpa(0.5 * c(s_station))); %wingbox width (m)
b_box = double(vpa(0.0942 * c(s_station))); %wingbox average height (m)
N = bm_new./(c_box.*b_box); %compressive stress in N/m


%Initial buckling
t_init = ((N .* b_box.^2)./(3.62 * E)).^(1/3) * 1000; %skin thickness in mm
cr_stress = N./(1000*t_init); %critical buckling stress of panel in N/mm^2 (sigma)
area_init = c_box .* t_init .* 1000; %cross-sectional area of skin with in mm^2


 Design panel size and stringer size
n_st = 16;
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

%Use initial buckling of panel to find skin thickness
t = 1000 .* ((N .* (b/1000).^2) ./ (3.62 * E)) .^ (1/3); %skin thickness with stringers added in mm
cr_s0 = N./(1000*t); %critical buckling stress in N/mm^2 (sigma0)
t_e = t + As./b; %panel effective thickness in mm (smeared thickness?)
As_e = b .* n .* t_e; %panel effective area in mm^2

%R_t = ts./t; %ratio of ts/t
R_t = ts./t(1)
%R_b = As./b./t; %ratio of As/bt
R_b = As./(b(1).*t(1))

% for readability, no data for any lower/higher than these inequalities 
i1 = find(R_t > 1.4);
i2 = find(R_t < 0.4);

R_t(i1) = 1.4;
R_t(i2) = 0.4;


%Getting sigma0/sigma ratio from Catchpole diagram
R_sigma = catchpole(R_b, R_t);

% i3 = find(isnan(R_sigma));
% for i = 1:length(i3)
%     if i3(i) == 1
%         R_sigma(i3(i)) = R_sigma(2);
%     else
%         R_sigma(i3(i)) = R_sigma(i3(i) - 1);
%     end
% end

cr_stress1 = R_sigma .* cr_s0;

%Getting FARRAR efficiency for each span station
F = farrarval(R_b, R_t)
mean_stress = F.*sqrt((N .* E*10^9)./ L) * 10^-6;
L = N(1) * E / (cr_stress1(1) * 10^6/F)^2 % length in m

%Panel weight calculations
rho = 2850; %density of AL2024 in kg/m^3

W = [W, As_e(1) * rho/10^6]; %weight per unit length in grams
end

%Plotting for skin thickness and stress distribution
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
ylim([0 1.7])
hold off

WEIGHT_SKIN = sum(As_e.* (b_h)/12)/1000000 * rho * 2

%Interaction between shear and compression
yspan = linspace(0,b_h, 1000);

shear_flow = load("shear_flow.mat"); %need to get this file from excel dont have yet
q0 = shear_flow.q0
q0_new = spline(yspan, q0, s_station) ./ 1000;
s_q0 = q0_new ./ t
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



Rib Design
% iterating through number of ribs for the weight optimisation
%n_ribs = 15:1:30;
% rib_array = 0.05:0.05:0.7;
rib_array = 0.44;
W_skin = ones(length(rib_array), 1)';
W_rib = ones(length(rib_array), 1)';
W_total = ones(length(rib_array), 1)';

for i = 1 : length(rib_array)
L = rib_array(i);
% no. of ribs
n_ribs = b_h/L;

% defining span stations
s_station = linspace(0,b_h, n_ribs);

%wing chord distribution starting at the wing-fuselage intersection
syms y
c(y)  = 2.85 - 0.3021*(y);
chord = double(vpa(c(s_station)));

% defining new bm matching the span stations
bm_new = spline(yspan, bm, s_station);

%solve this problem of bm_new at the end going to the negative which causes
%problems of having complex number later on

c_box = double(vpa(0.5 .* chord)); %wingbox width
b_box = double(vpa(0.0942 .* chord)); %wingbox average height (m)
N = bm_new./(c_box.*b_box); %compressive stress in N/m

%Rib thickness
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

rho = 2850; %density of skin material in kg/m^3
rho_r = 2570; %density of rib material in kg/m^3

W_rib(i) = sum((1/3) * rho * ( (1/(4*FARRAR^2))^(1/3) + (2/FARRAR^2)^(1/3)) .* ( (b_box .* t_rib_b./1000 .* N ) ./ E ).^(1/3)); %sum(t_e) * rho/10^6;
W_skin(i) = mean((rho/FARRAR) .* sqrt((N .* L) ./ E))*9.81;
W_total(i) = W_skin(i) + W_rib(i);

end

WEIGHT_RIB = sum(ones(1,12)./ 1000 .* (c_box .* b_box)) .* rho_r

%Plotting
% %plotting weight opt
% figure
% hold on
% plot(rib_array, W_skin, '-b', 'LineWidth', 1.5)
% plot(rib_array, W_rib, '-r', 'LineWidth', 1.5)
% plot(rib_array, W_total, '-k', 'LineWidth', 1.5)
% plot(0.45, W_total(9), 'xg', 'LineWidth', 1.5)
% hold off
% grid minor
% xlabel("Rib spacing (m)")
% ylabel("Weight per unit skin surface area $(N/m^2)$")
% ylabel(get(get(gca,'YLabel'),'String'),'Interpreter','latex');
% lg = legend('Weight of Skin + Stringers', 'Weight of Ribs', 'Total weight', 'Chosen rib spacing');
% % Set global properties for all axes
% set(findobj(gcf, 'type', 'axes'),'FontSize', 13, 'LineWidth', 1);
% set(findobj(gcf, 'type', 'line'), 'LineWidth', 1.5);
% set(lg, 'Interpreter','latex');
% xlabel(get(get(gca,'XLabel'),'String'),'Interpreter','latex');
% box on
% 
% 
% %plotting stress variation of skin-stringer panels with rib spacing
% figure
% hold on
% plot(rib_array, sigma_crit, '-b', 'LineWidth', 1.5)
% plot(linspace(0,max(rib_array),length(rib_array)),sigma_y*10^6.*ones(length(rib_array),1),'k--', 'LineWidth', 1.5)
% hold off
% grid minor
% xlabel("Rib spacing (m)")
% ylabel("Maximum stress exerted on the skin-stringer panel (N/m^2)", Interpreter="latex")
% lg = legend('Maximum stress', 'Yield stress of skin-stringer material');
% % Set global properties for all axes
% set(findobj(gcf, 'type', 'axes'),'FontSize', 13, 'FontWeight', 'Bold', 'LineWidth', 1);
% set(findobj(gcf, 'type', 'line'), 'LineWidth', 1.5);
% set(lg, 'Interpreter','latex');
% xlabel(get(get(gca,'XLabel'),'String'),'Interpreter','latex');
% ylabel(get(get(gca,'YLabel'),'String'),'Interpreter','latex');


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



D-section design
N = 12; % number of pseudo-ribs;
t_guess = 1 ;% mm
t1 = t_guess.*ones(1,N); % guess d cell thickness array , mm

Shear calculation
semi_span =  b_h; % wing semi span outside fuselage, m
bref = 0:(semi_span/N):semi_span - (semi_span/N);
E =  79.8*10^3; %Young's modulus, MPa

box_h = b_box(1) ; % wing box height, m
tresca_stress = 320 % yield shear stress for Al 2090, MPa

a = semi_span/(N+1) % d cell span length, m

%wing chord distribution starting at the wing-fuselage intersection
syms y
c(y) = 2.85 - 0.3021*(y);

b = double(0.15 .* c(bref)); % width of d cell, m

%b = double(0.5.*(pi.*3.*((box_h/2) + d_cell_width) - sqrt((3.*d_cell_width + (box_h/2).*(d_cell_width + 3.*(box_h/2)))))) % curved panel perimeters array , m - assume it's an ellipse, can't be bothered to get exact aerofoil perimeter

R = double(0.5.*((b.^2)./(box_h/2) + (box_h/2)^2./b)); % average radius of curvatures array, m

a/b(1)
b(1)/a

coeff_1 = b./sqrt(R(1).*(t1./1000))  % for a > b
coeff_2 = a./sqrt(R(1).*(t1./1000)) % for a < b


Getting Ks from plot
x = [0, 2, 4, 6 ,8 ,10, 12, 14];
y = [9, 12, 17, 23, 29, 37, 45, 54];

Ks = abs(spline(x,y,coeff_1))


tau_crit = Ks(1).*E.*(t1(1)./(b(1).*1000)).^2 % critical buckling stress array, MPa, we want tau_crit = tresca_stress

t_req = sqrt(tresca_stress./(Ks.*E)).*b*1000


Plotting
figure
set(findobj(gcf,'type','axes'),'FontSize',20,'FontWeight','Bold', 'LineWidth', 1.5);
plot(bref, t_req, 'xk', LineWidth=1.5)
grid minor
xlabel('Span length (m)', Interpreter="latex", FontSize=14)
ylabel("Pseudo-ribs thickness (mm)", Interpreter="latex", FontSize=12)
xlim([0 5.5])
box on
WEIGHT_DSECCTION = sum(t_req./1000 .* (b .* box_h .* (2/3)) .* 2570)

WING PLOT
semispan_vt = 11.32/2;
n = 500;
rib_no = 14;
stringer_no = 17;

lambda = 0.4;
sweep_LE = 40;
wing_root_vt = 2.85;
wing_tip_vt = 1.14;

fspar_pos = 0.2;
rspar_pos = 0.7;

rib_pos = linspace(0,semispan_vt,rib_no);

yspan_vt = linspace(0,semispan_vt,n);

LE_x = yspan_vt./tand(90-sweep_LE);

root_x = linspace(0,wing_root_vt,n);
root_y = zeros([1,n]);

tip_x = (semispan_vt/tand(90-sweep_LE))+linspace(0,wing_tip_vt,n);
tip_y = semispan_vt.*ones([1,n]);

TE_x = linspace(wing_root_vt,max(tip_x),n);

F_spar_x = linspace(fspar_pos*wing_root_vt,(semispan_vt/tand(90-sweep_LE))+fspar_pos*wing_tip_vt,n);
R_spar_x = linspace(rspar_pos*wing_root_vt,(semispan_vt/tand(90-sweep_LE))+rspar_pos*wing_tip_vt,n);

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
    line([rib_pos(i)/(tand(90-sweep_LE)),rib_pos(i)/(tand(90-sweep_LE))+rspar_pos*(wing_root_vt+rib_pos(i)*(wing_root_vt*(lambda-1))/(semispan_vt))],[rib_pos(i),rib_pos(i)],'Color','#0000FF','LineWidth',1.8)
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

ylim([0,7])
set(findobj(gcf, 'type', 'axes'),'FontSize', 13, 'FontWeight', 'Bold', 'LineWidth', 0.05);
set(findobj(gcf, 'type', 'line'), 'LineWidth', 1.5);
xlabel(get(get(gca,'XLabel'),'String'),'Interpreter','latex');
ylabel(get(get(gca,'YLabel'),'String'),'Interpreter','latex');
lg = legend;
set(lg, 'Interpreter','latex');
ylim([0,6])
axis equal
grid on
xlabel('X Position (m)')
ylabel('Y Position (m)')
hold off
