clear
clc

bluez=[0 0 1];
redz=[1 0 0];
greenz=[0 0.75 0];
lilaz=[0.75 0 0.75];
orangz=[1 0.5 0];
colorz={bluez, redz, greenz, orangz, lilaz};
density=2780;
Vtail_W = 688; 
n=450; %Number of discretizations
cbar=3.3267; 
lambda=0.7; 
ctip=2.7115;
croot=ctip/lambda;
l_unswept=4.99; 
sweep=35; 
fs_c=0.2; %Front spar chord percentage
rs_c=0.7; %Rear spar chord percentage
LF=3.75; 

%Solving moment balance for OEI case to find lift force
thrust=81737; 
engine_y= 2.0356/2+1.02; 
ac_xcg=18.80; 
ac_xcg_alt =19.2; 
vert_le_x=29.71; 

L=(engine_y*thrust)/(cbar/4+vert_le_x- ac_xcg);
L_alt=(engine_y*thrust)/(cbar/4+vert_le_x- ac_xcg_alt);


l=l_unswept/cosd(sweep);
dl=flip(linspace(0,l,n+1)); 

L0=(4*L)/(pi*l);
dL=L0*sqrt(1- (dl/l).^2); % Sectional lift profile

L0_alt=(4*L_alt)/(pi*l);
dL_alt=L0_alt*sqrt(1- (dl/l).^2);

figure
box on
set(findobj(gcf,'type','axes'),'FontName','Palatino','FontSize',12,'FontWeight','Bold', 'LineWidth', 1.5);
hold on
grid on
plot(dl,dL,'-','Color',colorz{1},'LineWidth',1.5);
plot(dl,dL_alt,'-','Color',colorz{2},'LineWidth',1.5);

xlabel('Position [m]')
ylabel('Sectional Lift [N/m]')
legend('MTOW Case', 'MZFW Case')
lift=trapz(dL,dl); 
hold off

W=LF*Vtail_W*9.81; 
yspan = linspace(0,l_unswept,100);
semispan = l_unswept;
S1 = 0.705;
S2 = S1*(lambda)^2;
WingBox_Volume = 1/3*(S1+S2+(S1*S2)^(1/2))*semispan;

Area_Span = (S1*(lambda^2-2*lambda+1))*(yspan./semispan).^2 + (S1*(-2+2*lambda)).*(yspan./semispan) + S1;
WingW_Span = ((W)/WingBox_Volume).*Area_Span;

Lsec=[0,(dL(1:end-1)+dL(2:end)).*(dl(1:end-1)-dl(2:end))/2];
SF=cumsum(Lsec); 
dM=[0,(SF(1:end-1)+SF(2:end)).*(dl(1:end-1)- dl(2:end))/2];
BM=cumsum(dM); 

Lsec_alt= [0,(dL_alt(1:end-1) + dL_alt(2:end)).*(dl(1:end-1)-dl(2:end))/2];
SF_alt=cumsum(Lsec_alt); 
dM_alt=[0,(SF_alt(1:end-1)+ SF_alt(2:end)).*(dl(1:end-1)-dl(2:end))/2];
BM_alt=cumsum(dM_alt); 

figure
box on
set(findobj(gcf,'type','axes'),'FontName','Palatino','FontSize',12,'FontWeight','Bold', 'LineWidth', 1.5);
hold on
grid on
plot(dl,SF,'-','Color',colorz{1},'LineWidth',1.5);
plot(dl,SF_alt,'-','Color',colorz{2},'LineWidth',1.5);
xlabel('Position [m]')
ylabel('Shear Force [N]')
legend('MTOW Case', 'MZFW Case');

figure
box on
set(findobj(gcf,'type','axes'),'FontName','Palatino','FontSize',12,'FontWeight','Bold', 'LineWidth', 1.5);
hold on
grid on
plot(dl,BM,'-','Color',colorz{1},'LineWidth',1.5);
plot(dl,BM_alt,'-','Color',colorz{2},'LineWidth',1.5);
xlabel('Position [m]')
ylabel('Bending Moment [Nm]')
legend('MTOW Case', 'MZFW Case');

save("bm_vtail",'BM'); %BENDING MOMENT OF VTAIL
save("sf_vtail",'SF'); %SHEAR FORCE OF VTAILL

c_swept=flip(linspace(croot,ctip,n+1))*cosd(sweep); 
flex_c=(fs_c+rs_c)/2;
dT=Lsec.*(flex_c-0.25).*c_swept; %Sectional torque (Nm)
T=cumsum(dT);

dT_alt=Lsec_alt.*(flex_c-0.25).*c_swept; %Sectional torque (Nm)
T_alt=cumsum(dT_alt);

figure
box on
set(findobj(gcf,'type','axes'),'FontName','Palatino','FontSize',12,'FontWeight','Bold', 'LineWidth', 1.5);
hold on
grid on
plot(dl,T,'-','Color',colorz{1},'LineWidth',1.5);
plot(dl,T_alt,'-','Color',colorz{2},'LineWidth',1.5);
xlabel('Position [m]')
ylabel('Torque [Nm]')
legend('MTOW Case', 'MZFW Case');

save("torque_vtail",'T'); 
maxvalue = [max(abs([SF,SF_alt])),max(abs([BM,BM_alt])),max(abs([T,T_alt]))]

clear; clc; close all;

%% Wing Geometry for Vertical Stabilizer
l_unswept = 6.8;  
sweep_angle = 34.1;  
l_swept = l_unswept / cosd(sweep_angle);  

yspan_data = load("yspan_vtail.mat"); 
yspan = yspan_data.yspan_v;
bm_data = load("bm_vtail.mat");  
bm_original = flip(bm_data.BM);  
sf_data = load("sf_vtail.mat");  
sf_original = flip(sf_data.SF); 
torque_data = load("torque_vtail.mat");  
torque_original = torque_data.T;  


span_stations = 0:0.64:l_swept;

%% Material Properties for Composite Laminate

sigma_xx_tensile = 2132;  
sigma_xx_comp = 1548;  
sigma_xy = 50;  

% Fiber and matrix properties (MPa)
sigma_f_tensile = 3500;  
sigma_m_tensile = 80;  
sigma_f_comp = 2500;  
sigma_m_comp = 120;  
tau_f = 60;  
tau_m = 50;  
E_f = 230000;  
E_m = 3000;  
vf = 0.60;  
vm = 0.40;  
sigma_t_1 = sigma_f_tensile * vf + sigma_m_tensile * vm;  
sigma_t_2 = 50;  
sigma_c_1 = sigma_f_comp * vf + sigma_m_comp * vm;  
sigma_c_2 = sigma_m_comp * (1 - 2 * sqrt(vf / pi));  
tau_12 = tau_m;  
E_1 = E_f * vf + E_m * vm;  % Longitudinal modulus (MPa)
E_2 = (vf / E_f + vm / E_m)^-1;  % Transverse modulus (MPa)

t_ply = 0.2; 

SF = spline(yspan, sf_original, span_stations);  % Shear force interpolation
BM = spline(yspan, bm_original, span_stations);  % Bending moment interpolation

syms x
c(x) = -0.6618 * x + 7.5;  % Chord function

% Compute chord values at defined spanwise stations
c_dis = double(vpa(c(span_stations)));

h_web = 0.0942 .* c_dis;  % Wingbox height (m)
c_web = 0.5 .* c_dis;  % Wingbox width (m)

frontspar = 0.2;
rearspar = 0.7;

% Flange dimensions
flange_web_ratio = 0.3;  
b_flange = flange_web_ratio * h_web;  

% Calculate spar moment of inertia (m^4)
Ixx_spar = (h_web(1).^3 .* (1.25 / 1000)) / 12 + ...
           2 * (((1.25 / 1000)^3 .* b_flange(1)) / 12 + ((1.25 / 1000) .* b_flange(1) .* h_web(1).^2) / 4);

%% Spar Cap Layer Calculations
n_sc_tensile = BM ./ ((sigma_xx_tensile * 1e6) * (t_ply / 1000) .* h_web .* b_flange); 
n_sc_comp = BM ./ ((sigma_xx_comp * 1e6) * (t_ply / 1000) .* h_web .* b_flange);  
n_sc_manu = ceil(n_sc_comp);  

% Manually adjusted number of plies for critical areas based on design criteria
n_sc_manu_new = [12, 12, 12, 12, 10, 10, 10, 10, 8, 8, 8, 8, 8];  

%% Critical Buckling Load (Nx_crit) Calculation

D11 = [446436, 338968, 250237];  
D22 = [113602, 81709.4, 56042.1];  
D12 = [59536.5, 44774.3, 32837.6];  
D66 = [66708.7, 50162.8, 36765.9];  
span_web = [0, 2, 4];
b_web = spline(span_stations, h_web, span_web);  
n_stiff = [1, 1, 1];  

b_des = b_web ./ n_stiff;  

N_x_crit = ((2 * pi^2) ./ b_des.^2) .* ...
           ((sqrt(D11 * 1e-3 .* D22 * 1e-3) + D12 * 1e-3 + 2 * D66 * 1e-3));  

N_x = BM ./ (h_web .* c_web);  % Applied normal force per unit length (N/m)

figure 
hold on
stairs([span_web 6], [N_x_crit N_x_crit(end)])
plot(span_stations, N_x)
legend('N_{x,crit}', 'N_{xx}')
ylabel('Buckling load (N)')
xlabel('Distance along span (m)')
xlim([0,6.8])
grid minor

n_sw = (0.5 .* SF)./((sigma_xy*10^6).*(t_ply/1000).*h_web)
n_sw_manu = ceil(n_sw)

D11 = [6211.04   1442.71];
D22 = [2707.51   942.208];
D12 = [2428.77   848.444];
D66 = [2574.27   891.553];

span_web = [0 1]

b_web = spline(span_stations, h_web, span_web);

n_stiff = [6 10]; 

b_des = b_web./n_stiff;

N_x_crit = ((2*pi^2)./b_des.^2).*((sqrt(D11.*10^-3.*D22.*10^-3) + D12.*10^-3 + 2.*D66.*10^-3)); % N/m

B = (b_web./0.5).*((D11.*10^-3)./(D22.*10^-3)).^(1/4)
phi = (sqrt(D11.*10^-3.*D22.*10^-3))./(D12.*10^-3 + 2.*D66.*10^-3)

ks = [12 11];     % from graph

N_xy_crit = ((ks.*pi^2)./b_des.^2).*(D11.*10^-3.*(D22.*10^-3).^3).^(1/4);

b = abs((frontspar-rearspar)*c_dis);
N_x_old = BM./(b.*h_web);
N_x_applied = spline(span_stations, N_x_old, span_web);


N_xy_old = SF./h_web; % N/m
N_xy_applied = spline(span_stations, N_xy_old, span_web);

web_buckle = (N_xy_applied./N_xy_crit).^2 + (N_x_applied./N_x_crit)

% this is the stress in the VTail in the xy direction
N_xy = SF./h_web; 

D11 = 713593;
D22 = 221598;
D12 = 92887.6;
D66 = 104726;

n_stringers = 4;
b = abs((frontspar-rearspar)*c_dis);
b_des = (abs((frontspar-rearspar)*c_dis))/n_stringers;
N_x_crit = ((2*pi^2)./b_des(1)^2).*((sqrt(D11*10^-3*D22*10^-3) + D12*10^-3 + 2*D66*10^-3)); % N/m
N_x_applied = BM(1)/(b(1)*h_web(1)) % N/m

skin_buckle = N_x_applied./N_x_crit

t_e = 3.25; % skin composite thickness in mm
no_ribs = 5;
l=6.8;
ribspace = l/no_ribs % rib spacing
span_rib = linspace(0,l,no_ribs) 

cdis_new = spline(span_stations, c_dis, span_rib);
b = abs((frontspar-rearspar)*c_dis);
b_new = spline(span_stations, b, span_rib);
h_web = 0.7*0.12.*cdis_new;
h_c = h_web;
Is = ((cdis_new.*(t_e./1000).^3)/12) + (cdis_new.*(t_e./1000).*(h_c./2).^2);
E = 133; 
BM_new = spline(span_stations, BM, span_rib);

% Calculate the crushing load
F_c = ((BM_new.^2).*ribspace.*h_c.*(t_e/1000).*cdis_new)./(2*E*10^9.*(Is.^2)); %(N)

N_x_applied = F_c./b_new

D11  = 203002;  
D22 = 50848.9; 
D12 = 3692.41; 
D66 = 6451.36; 

b_des = b_new;
N_x_crit = ((2*pi^2)./b_des.^2).*((sqrt(D11*10^-3*D22*10^-3) + D12*10^-3 + 2*D66*10^-3)); % N/m

rib_buckle = N_x_applied./N_x_crit

figure
set(findobj(gcf,'type','axes'),'FontName','Palatino','FontSize',12,'FontWeight','Bold', 'LineWidth', 1.5);
plot(span_stations, n_sc_comp, '--r', 'LineWidth', 1.5)
hold on
stairs(span_stations ,n_sc_manu, 'b', 'LineWidth', 1.5)
stairs(span_stations,n_sc_manu_new, 'k','LineWidth', 1.5)
grid on
box on
legend('Theoretical number ${n_0}$', 'Step theoretical number ${n_0}$', 'Final step number ${n_0}$ ','FontSize',12,'FontWeight','Bold','interpreter','latex')
xlabel('Distance along span [m]','FontSize',14,'interpreter','latex')
ylabel('Number of ${n_0}$ in the spar caps','FontSize',14,'interpreter','latex')
ylim([0 15])
xlim([0,6.8])
grid minor
box on

figure
set(findobj(gcf,'type','axes'),'FontName','Palatino','FontSize',12,'FontWeight','Bold', 'LineWidth', 1.5);
plot(span_stations, n_sw, '--r', 'LineWidth', 1.5)
hold on
stairs(span_stations ,n_sw_manu, 'b', 'LineWidth', 1.5)
grid minor
box on
xlim([0,6.8]);
legend('Theoretical number ${n_{45}}$', 'Manufactured number ${n_{45}}$','FontName','Palatino', 'FontSize',12,'interpreter','latex')
xlabel('Distance along span [m]','FontSize',14,'interpreter','latex')
ylabel('Number of ${n_{45}}$ in the spar caps','FontSize',14,'interpreter','latex')

syms x

ht.b = 6.8*2;
AR_v = 1.6;
sweep_v_quarter = 34.1;
lambda_v = 0.4;

sweep_v_LE = atand(tand(sweep_v_quarter) - (4/AR_v)*(((0-25)/100)*((1 - lambda_v)/(1 + lambda_v))));
sweep_v_TE = 14;

% Rib Spacing
rib_spacing = 1.6424; 
num_ribs = floor(6.8 / rib_spacing); 

% Generate rib positions from root to tip
L_output = linspace(0, num_ribs * rib_spacing, num_ribs+1); 

%Line representing LE
LE(x) = 6.8 - (tand(sweep_v_LE))*x;

%Line representing TE
TE(x) = -(tand(sweep_v_TE))*x;

figure

xlabel('Spanwise length (m)','FontSize',13,'interpreter','latex')
ylabel('Chordwise length (m)','FontSize',13,'interpreter','latex')
grid on
box on
hold on
ax = gca;
ax.FontSize = 15;
set(gcf,'units','inches','position',[1,1,8,6])

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% --- Plotting the actual wing ---
fplot(LE(x) , [0,ht.b/2],'k','LineWidth', 1.5) 
fplot(TE(x) , [0,ht.b/2],'k','LineWidth', 1.5) 
plot([ht.b/2 ht.b/2], [TE(ht.b/2) LE(ht.b/2)],'k','LineWidth', 1.5) 

h_spars(1) = fplot(TE(x)+0.35*(LE(x)-TE(x)), [0,ht.b/2], 'r','LineWidth', 1.5); 
h_spars(2) = fplot(LE(x)-0.1*(LE(x)-TE(x)), [0,ht.b/2], 'r','LineWidth', 1.5); 

b = 230;
for i = 1:16
    y1 = LE(x)-0.1*(LE(x)-TE(x));
    y2 = TE(x)+0.35*(LE(x)-TE(x)) + i*(b/1000);
    eqn = y1 == y2;
    f(i) = double(vpa(solve(eqn, x)));

    if f(i) > ht.b/2
        f(i) = ht.b/2;
    end

    h_stringers(i) = fplot(TE(x)+0.35*(LE(x)-TE(x)) + i*(b/1000), [0, f(i)], 'g','LineWidth', 0.8);

    y_start = subs((TE(x)+0.35*(LE(x)-TE(x)) + i*(b/1000)), 0);
    y_end = subs((TE(x)+0.35*(LE(x)-TE(x)) + i*(b/1000)), f(i));

    length_stringers(i) = double(vpa(sqrt( (f(i)-0)^2 + (y_end - y_start)^2) , 2 ));
end


for i = 1:length(L_output)
    sum1 = L_output(i);

    y = linspace(LE(sum1) - 0.1*(LE(sum1) - TE(sum1)), TE(sum1) + 0.35*(LE(sum1) - TE(sum1)), 100);
    x = sum1 * ones(1, length(y)); 

    % Plot the rib
    h_ribs(i) = plot(x, y, 'b', 'LineWidth', 1.5);
    length_ribs(i) = double(vpa(sqrt((x(100)-x(1))^2 + (y(100) - y(1))^2), 2));
end
legend([h_ribs(1), h_spars(1), h_stringers(1)], 'VT-Ribs', 'VT-Spars', 'VT-Stringers', 'Location', 'best')

hold off
