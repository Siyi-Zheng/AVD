%VTail Composite Design

% housekeeping

clear,clc, close all


% wing geometry
l_unswept=6.8; %Length of vertical stabiliser (m) without sweep (aka in the z direction)
sweep=34.1; %Stabiliser sweep angle (deg)
l=l_unswept/cosd(sweep); %Swept length of stabiliser
% dl=flip(linspace(0,l,499)); %Discretizing wing into 500 sections AND flipping it so it goes from tip to root

yspan = load("yspan_vtail.mat"); %array of y coords spanwise
yspan = yspan.yspan_v;

bm_original = load("bm_vtail.mat")
bm_original = flip(bm_original.BM)
sf_original = load("sf_vtail.mat")
sf_original = flip(sf_original.SF)
torque_original = load("torque_vtail.mat")
torque_original = torque_original.T

spanstation = [0:0.64:l]

sigma_xx_tensile = 2132; % MPa 
sigma_xx_comp = 1548; % MPa
sigma_xy = 50; % MPa

sigma_f_tensile = 3500; % MPa
sigma_m_tensile = 80; % MPa
sigma_f_comp = 2500; % MPa
sigma_m_comp = 120; % MPa
tau_f = 60; % MPa
tau_m = 50; % MPa
E_f = 230000; % MPa
E_m = 3000; % MPa

vf = 0.60;
vm = 0.40;

sigma_t_1 = sigma_f_tensile*vf + sigma_m_tensile*vm; % 2132 MPa
sigma_t_2 = 50; % MPa
sigma_c_1 = sigma_f_comp*vf + sigma_m_comp*vm; % 1548 MPa
sigma_c_2 = sigma_m_comp*(1 - 2*sqrt(vf/pi)); % 14.76 MPa
tau_12 = tau_m; % 50 MPa
E_1 = E_f*vf + E_m*vm; % 139200 MPa
E_2 = (vf/E_f + vm/E_m)^-1; % 7500 MPa

t_ply = 0.2; % mm, ply thickness

%discritising bm and sf into spanstations
SF = spline(yspan, sf_original, spanstation);
BM = spline(yspan, bm_original, spanstation);

syms x
c(x) = -0.6618*x + 7.5;

for i = 1:length(spanstation)
    cdis(i) = double(vpa(c(spanstation(i))));
end

h_web = 0.0942 .* cdis; % wingbox height
c_web = 0.5 .* cdis; % wingbox width

frontspar = 0.2;
rearspar = 0.7;

flange_web_ratio = 0.3; % GUESS
b_flange = flange_web_ratio*h_web;
Ixx_spar = (h_web(1).^3.*(1.25/1000))/12 + 2*(((1.25/1000)^3.*b_flange(1))/12 + ((1.25/1000).*b_flange(1).*h_web(1).^2)/4); % m^4


Spar cab

n_sc_tensile = BM./((sigma_xx_tensile*10^6).*(t_ply/1000).*h_web.*b_flange);
n_sc_comp = BM./((sigma_xx_comp*10^6).*(t_ply/1000).*h_web.*b_flange);
n_sc_manu = ceil(n_sc_comp);
% we found that these layers are not enough to satisfy the Nx_crit so we
% will manually increase the amount of plier number ourself

n_sc_manu_new = [12 12 12 12 10 10 10 10 8 8 8 8 8];


% checking for Nx_crit
% D11 = [446436   338968   250237]
% D22 = [113602   81709.4   56042.1]
% D12 = [59536.5   44774.3   32837.6]
% D66 = [66708.7   50162.8   36765.9]
D11 = 446436.0 
D22 = 113602.0 
D12 = 59536.5 
D66 = 66708.7 

span_web = [0 2 4];

b_web = spline(spanstation, h_web, span_web); % in meters
%adding stiffeners to get it to not buckle. Apparently it buckles all along
%span
n_stiff = [1 1 1]; 

b_des = b_web./n_stiff;

% this is the Nx_crit stress
N_x_crit = ((2*pi^2)./b_des.^2).*((sqrt(D11.*10^-3.*D22.*10^-3) + D12.*10^-3 + 2.*D66.*10^-3)) % N/m

% this is the stress in the VTail in x direction
N_x = BM ./ (h_web .* c_web)

Plotting to see which spanstation exceeds the Nx_crit stress
figure 
hold on
stairs([span_web 6], [N_x_crit N_x_crit(end)])
plot(spanstation, N_x)
legend('Nx_crit', 'N')
ylabel('Buckling load (N)')
xlabel('Distance along span (m)')
xlim([0,6.8])
grid minor

Spar Web
n_sw = (0.5 .* SF)./((sigma_xy*10^6).*(t_ply/1000).*h_web)
%n_sw_manu = [ceil(n_sw(1:11)) floor(n_sw(12:13))];
n_sw_manu = ceil(n_sw)

D11 = 6211.04;
D22 = 2707.51;
D12 = 2428.77;
D66 = 2574.27;

span_web = [0 1]

b_web = spline(spanstation, h_web, span_web);

n_stiff = [6 10]; 

b_des = b_web./n_stiff;

N_x_crit = ((2*pi^2)./b_des.^2).*((sqrt(D11.*10^-3.*D22.*10^-3) + D12.*10^-3 + 2.*D66.*10^-3)); % N/m

B = (b_web./0.5).*((D11.*10^-3)./(D22.*10^-3)).^(1/4)
phi = (sqrt(D11.*10^-3.*D22.*10^-3))./(D12.*10^-3 + 2.*D66.*10^-3)

ks = [12 11];     % from graph

N_xy_crit = ((ks.*pi^2)./b_des.^2).*(D11.*10^-3.*(D22.*10^-3).^3).^(1/4);

b = abs((frontspar-rearspar)*cdis);
N_x_old = BM./(b.*h_web);
N_x_applied = spline(spanstation, N_x_old, span_web);


N_xy_old = SF./h_web; % N/m
N_xy_applied = spline(spanstation, N_xy_old, span_web);


Check interaction
web_buckle = (N_xy_applied./N_xy_crit).^2 + (N_x_applied./N_x_crit)

% this is the stress in the VTail in the xy direction
N_xy = SF./h_web; % N/m


Composite Skin
D11 = 713593;
D22 = 221598;
D12 = 92887.6;
D66 = 104726;

n_stringers = 4;
b = abs((frontspar-rearspar)*cdis);
b_des = (abs((frontspar-rearspar)*cdis))/n_stringers;
N_x_crit = ((2*pi^2)./b_des(1)^2).*((sqrt(D11*10^-3*D22*10^-3) + D12*10^-3 + 2*D66*10^-3)); % N/m
N_x_applied = BM(1)/(b(1)*h_web(1)) % N/m

skin_buckle = N_x_applied./N_x_crit


Composite Rib
t_e = 3.25; % skin composite thickness in mm
no_ribs = 5;
ribspace = l/no_ribs%rib spacing
span_rib = linspace(0,l,no_ribs) 

cdis_new = spline(spanstation, cdis, span_rib);
b = abs((frontspar-rearspar)*cdis);
b_new = spline(spanstation, b, span_rib);
h_web = 0.7*0.12.*cdis_new;
h_c = h_web;
Is = ((cdis_new.*(t_e./1000).^3)/12) + (cdis_new.*(t_e./1000).*(h_c./2).^2);
E = 133; % GPa
BM_new = spline(spanstation, BM, span_rib);

% Calculate the crushing load
F_c = ((BM_new.^2).*ribspace.*h_c.*(t_e/1000).*cdis_new)./(2*E*10^9.*(Is.^2)); %(N)

N_x_applied = F_c./b_new

D11  = 203002;  % Nmm
D22 = 50848.9; % Nmm
D12 = 3692.41; % Nmm
D66 = 6451.36; % Nmm

b_des = b_new;
N_x_crit = ((2*pi^2)./b_des.^2).*((sqrt(D11*10^-3*D22*10^-3) + D12*10^-3 + 2*D66*10^-3)); % N/m

rib_buckle = N_x_applied./N_x_crit

Plotting for number of pliers
figure
set(findobj(gcf,'type','axes'),'FontName','Palatino','FontSize',12,'FontWeight','Bold', 'LineWidth', 1.5);
plot(spanstation, n_sc_comp, '--r', 'LineWidth', 1.5)
hold on
stairs(spanstation ,n_sc_manu, 'b', 'LineWidth', 1.5)
stairs(spanstation,n_sc_manu_new, 'k','LineWidth', 1.5)
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
plot(spanstation, n_sw, '--r', 'LineWidth', 1.5)
hold on
stairs(spanstation ,n_sw_manu, 'b', 'LineWidth', 1.5)
grid minor
box on
xlim([0,6.8]);
legend('Theoretical number ${n_{45}}$', 'Manufactured number ${n_{45}}$','FontName','Palatino', 'FontSize',12,'interpreter','latex')
xlabel('Distance along span [m]','FontSize',14,'interpreter','latex')
ylabel('Number of ${n_{45}}$ in the spar caps','FontSize',14,'interpreter','latex')


Plot of TAIL WING
syms x

ht.b = 6.8*2;
AR_v = 1.6;
sweep_v_quarter = 34.1;
lambda_v = 0.4;

sweep_v_LE = atand(tand(sweep_v_quarter) - (4/AR_v)*(((0-25)/100)*((1 - lambda_v)/(1 + lambda_v))));
sweep_v_TE = 14;

% Rib Spacing
rib_spacing = 1.6424; % Distance between consecutive ribs (meters)
num_ribs = floor(6.8 / rib_spacing); % Total number of ribs along span

% Generate rib positions from root to tip
L_output = linspace(0, num_ribs * rib_spacing, num_ribs+1); % Rib positions

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
grid on
box on
ax = gca;
ax.FontSize = 15;
set(gcf,'units','inches','position',[1,1,8,6])
% 
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

%plotting the actual wing
fplot(LE(x) , [0,ht.b/2],'k','LineWidth', 1.5) %plotting the leading edge
fplot(TE(x) , [0,ht.b/2],'k','LineWidth', 1.5) %plotting the trailing edge
plot([ht.b/2 ht.b/2], [TE(ht.b/2) LE(ht.b/2)],'k','LineWidth', 1.5) %plot the line to connect LE and TE

%plotting the spars
fplot(TE(x)+0.35*(LE(x)-TE(x)), [0,ht.b/2], 'r','LineWidth', 1.5) %near TE spars
fplot(LE(x)-0.1*(LE(x)-TE(x)), [0,ht.b/2], 'r','LineWidth', 1.5) %near LE spars

%plotting the stringers
b = 230;
for i = 1:11
    y1 = LE(x)-0.1*(LE(x)-TE(x));
    y2 = TE(x)+0.35*(LE(x)-TE(x)) + i*(b/1000);
    eqn = y1 == y2;
    f(i) = double(vpa(solve(eqn, x)));

    if f(i) > ht.b/2
        f(i) = ht.b/2;
    end

    fplot(TE(x)+0.35*(LE(x)-TE(x)) + i*(b/1000), [0, f(i)], 'g','LineWidth', 0.8)

    y_start = subs((TE(x)+0.35*(LE(x)-TE(x)) + i*(b/1000)), 0);
    y_end = subs((TE(x)+0.35*(LE(x)-TE(x)) + i*(b/1000)), f(i));

    length_stringers(i) = double(vpa(sqrt( (f(i)-0)^2 + (y_end - y_start)^2) , 2 ));

end

%plotting the ribs
for i = 1:length(L_output)
    sum1 = L_output(i); % Current rib location along the span

    % Define y-coordinates between the spars
    y = linspace(LE(sum1) - 0.1*(LE(sum1) - TE(sum1)), TE(sum1) + 0.35*(LE(sum1) - TE(sum1)), 100);
    x = sum1 * ones(1, length(y)); % Keep x constant at rib location

    % Plot the rib
    plot(x, y, 'b', 'LineWidth', 1.5);

    % Calculate the actual rib length
    length_ribs(i) = double(vpa(sqrt((x(100)-x(1))^2 + (y(100) - y(1))^2), 2));
end


Weight
% Weights

boxw_distb = b_new(1:5);
t_rb = 1.5;

stringers_weight = 2.*1550.*(120/10^6).*length_stringers;
stringers_weight_total = sum(stringers_weight); % KG

ribs_weight = 1550.*length_ribs.*boxw_distb.*(t_rb/1000);
ribs_weight_total = sum(ribs_weight); % KG

skin_weight_total = 2*1550*10.25*(2.5/1000); % KG

spar_weight = (1550*(2*(1.3125/1000)*0.16 + h_web(1)*(1.3125/1000)))*(6.11 + 7.03)

