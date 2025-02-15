%Dependencies: StringerCounterFunc.m and SweepCalcFunc.m
clear
clc
bluez=[0 0 1];
redz=[1 0 0];
greenz=[0 0.75 0];
lilaz=[0.75 0 0.75];
orangz=[1 0.5 0];
colorz={bluez, redz, greenz, orangz, lilaz};
density= 1565; %of tailplane material
Vtail_W = 1900; %Weight of tailplane kg


n=451; %Number of discretizations

cbar= 3.7; %Mean aerodynamic chord of vertical stabiliser (m)
lambda=0.55; %Taper ratio
ctip= 1.5; %tip chord length
croot=ctip/lambda;
l_unswept= 6.5; %Length of vertical stabiliser (m) without sweep (in the z direction)
sweep=34.1; %Stabiliser sweep angle (deg)

fs_c=0.2; %Front spar chord percentage
rs_c=0.7; %Rear spar chord percentage

LF=3.75; %load factor

%Solving moment balance for OEI case to find lift force
thrust=322000; %Engine thrust (N)
engine_y= 3.7; %Engine cg y coordinate (m)
ac_xcg=36; %Aircraft overall x_cg (m), using MTOW as aftmost
ac_xcg_alt =35.6; %x_cg (m) for ZW case
vert_le_x= 70.3; %x coordinate of vertical stabiliser leading edge at root (m)

L=(engine_y*thrust)/(cbar/4+vert_le_x- ac_xcg);
L_alt=(engine_y*thrust)/(cbar/4+vert_le_x- ac_xcg_alt);


%DOING DISTRIBUTION ALONG SWEPT COORDINATE
l=l_unswept/cosd(sweep); %Swept length of stabiliser
dl=flip(linspace(0,l,n+1)); %Discretizing wing into 1000 sections AND flipping it so it goes from tip to root

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
lift=trapz(dL,dl); % (might be slightly different due to integration errors)
hold off


%--------------------------Inertial loads------------------------------
%Not needed in vertical stabiliser case (perpendicular to lift)
%But I guess it's nice to have in case of structural calculations
W=LF*Vtail_W*9.81; %Vertical stabiliser weight (N)
% wSecMax=(2*W)/(l_unswept*(ctip/croot+ 1)); %Maximum sectional weight i.e. at the root (in N/m)
% wSecMin=(ctip/croot)*wSecMax; %Minimum sectional weight i.e. at the tip (in N/m)
% 
% %Weight distribution along horizontal stabiliser (straight line equation):
% wSec=((wSecMin-wSecMax)/l_unswept).*flip(dl)+wSecMax;
% 
% %Plotting weight distribution
% figure
% box on
% set(findobj(gcf,'type','axes'),'FontName','Palatino','FontSize',12,'FontWeight','Bold', 'LineWidth', 1.5);
% hold on
% grid on
% plot(flip(dl),wSec,'LineWidth',1.5)
% xlabel('Spanwise Position [m]')
% ylabel('Sectional Weight [N/m]')
% 
% hold off
%Weight of ONE Wing in N
yspan = linspace(0,l_unswept,100);
semispan = l_unswept;
S1 = 6.96;%Area wingbox  
S2 = S1*(lambda)^2;
WingBox_Volume = 1/3*(S1+S2+(S1*S2)^(1/2))*semispan;

Area_Span = (S1*(lambda^2-2*lambda+1))*(yspan./semispan).^2 + (S1*(-2+2*lambda)).*(yspan./semispan) + S1;
WingW_Span = ((W)/WingBox_Volume).*Area_Span;



% 
% plot(yspan,WingW_Span,'r--')
% set(findobj(gcf,'type','axes'),'FontName','Palatino','FontSize',12,'FontWeight','Bold', 'LineWidth', 1.5);
% legend('Sectional Weight')
% xlabel('Position [m]')
% ylabel('Sectional Weight [N/m]')
% title('Vertical Tail: Load Case 2')
% grid on




%-----------------Shear Force and Bending Moment-----------------------
Lsec=[0,(dL(1:end-1)+dL(2:end)).*(dl(1:end-1)-dl(2:end))/2];
SF=cumsum(Lsec); %Hehe funni
dM=[0,(SF(1:end-1)+SF(2:end)).*(dl(1:end-1)- dl(2:end))/2];
BM=cumsum(dM); %funni 2.0

Lsec_alt= [0,(dL_alt(1:end-1) + dL_alt(2:end)).*(dl(1:end-1)-dl(2:end))/2];
SF_alt=cumsum(Lsec_alt); %Hehe funni
dM_alt=[0,(SF_alt(1:end-1)+ SF_alt(2:end)).*(dl(1:end-1)-dl(2:end))/2];
BM_alt=cumsum(dM_alt); %funni 2.0

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


%-----------------------------Torque-----------------------------------
c_swept=flip(linspace(croot,ctip,n+1))*cosd(sweep); %Discretize Spanwise chord distribution into 1000 sections, from tip to root, swept length (used only for Torque calculation)
flex_c=(fs_c+rs_c)/2; %Chordwise percentage Position of flexural axis (assumed to be halfway between both spars)
%No pitching/M0 (symmetric airfoil) loading and no weight/inertial loading (vertical stabiliser)
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

save("torque_vtail",'T'); %effective panel thick

% Vertical stabiliser planform dimensions:
S_v = 29.1; % Reference area in m^2 (Change based on design)
AR_v = 1.6; % Aspect ratio of the vertical stabilizer
lambda_v = 0.55; % Taper ratio (tip chord / root chord)
sweep_v = 34.1; % Leading-edge sweep angle of vertical stabilizer (in degrees)

bperp_v = sqrt(AR_v * S_v); % Vertical span (height of vertical stabilizer)
b_v = bperp_v / cosd(sweep_v); % Actual slant height considering sweep
s_v = bperp_v ; % Semispan (half of the vertical stabilizer height)

% Spanwise points from root to tip
yspan_v = linspace(0, bperp_v, 452); 

% Spanwise position at root and tip
yroot_v = 0; % Root (base of the vertical stabilizer)
ytip_v = b_v ; % Tip (top of the vertical stabilizer)

% Save the spanwise coordinate data
save('yspan_vtail', 'yspan_v')

maxvalue = [max(abs([SF,SF_alt])),max(abs([BM,BM_alt])),max(abs([T,T_alt]))]
