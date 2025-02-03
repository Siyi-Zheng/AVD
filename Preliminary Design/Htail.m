clear
clc
bluez=[0 0 1];
redz=[1 0 0];
greenz=[0 0.75 0];
lilaz=[0.75 0 0.75];
orangz=[1 0.5 0];
colorz={bluez, redz, greenz, orangz, lilaz};

n=1000; % Number of discretisations you want
%% Getting geometries
S_ref = 482;

% Horizontal stabiliser planform dimensions:
S=58; % Reference area in m^2
AR=5.8; % Aspect ratio
lambda=0.4; % Taper ratio
bperp=sqrt(AR*S); % Span (perpendicular to fuselage)
sweep=35.155; % Flexural axis sweep angle of horizontal stabiliser (in degrees)
b=bperp/cosd(sweep); % Length of wing
s=bperp/2; % Semispan of wing
yspan = linspace(0,bperp/2,1000);% Getting spanwise points from root to tip
yroot=0; % Spanwise position at root
ytip=b/2; % Spanwise position at tip (=semispan)
save('yspan_htail', 'yspan')
% Getting distances of wing+h.stab ac and cg of aircraft (in m):
xcg=36;
xcg_zf = 35.6;
xw=35.3;
% Getting chord distribution:
croot=(2*S)/(bperp*(1+lambda)); % Chord at root (in m)
ctip=lambda*croot; % Tip chord in m
c=((ctip-croot)/s).*yspan + croot;
xh=56; % Assuming ac is at 1/4 chord point


%% Getting total lift of horizontal stabilisers

% Idea is that lift from h.stab needs to balance moment from lift of wing
% about CG of aircraft - assume that all lift to balance weight comes from
% the wing and the lift from the h.stab only exists when considering
% moments

% Getting lift from wings:
g=9.81; % Gravitational acceleration
ULF=3.75; % Ultimate load factor (n=1.5*2.5=3.75)
MTOW=353000*g; % MTOW in N
MZFW = 162976*g; %MZFW
Lw=ULF*MTOW; % Lift from wings in N
Lw_zf = ULF*MZFW;

% Getting moment (lever) arms:
leverw=xcg-xw; % Moment arm of wing
leverh=xh-xcg; % Moment arm of horizontal stabiliser

leverw_zf=xcg_zf-xw; % Moment arm of wing
leverh_zf=xh-xcg_zf; % Moment arm of horizontal stabiliser

VD=214.06; % Dive speed (m/s, IAS)
VA = 171.25;
Cm0=-0.131; % Zero lift pitching moment coefficient
M0w=0.5*1.225*VD^2*S_ref*Cm0; % Zero lift pitching moment from wings
M0w_A=0.5*1.225*VA^2*S_ref*Cm0;
% Using moment balance about aircraft CG to find lift from h.stab:
Lh=(Lw*leverw-M0w)/leverh; % Lift of horizontal stabiliser in N
Lh_zf=(Lw_zf*leverw_zf-M0w)/leverh_zf; % Lift of horizontal stabiliser in N

%% Getting lift distribution of horizontal stabiliser

% SINCE SYMMETRICAL, ONLY NEED TO CONSIDER ONE WING:
Lh2=0.5*Lh; % Lift of one horizontal stabiliser in N
Lh2_zf=0.5*Lh_zf;
%s=b/2; % Semispan of horizontal stabilisers

% Assuming elliptical distribution: total lift of one horizontal stabiliser
% is 1/4 of an ellipse so use 1/4 of area of an ellipse and the semi-span
% to find sectional lift at the root:
LSecMax=(4*Lh2)/(pi*s); % in N
LSecMax_zf=(4*Lh2_zf)/(pi*s);

LSec=sqrt((LSecMax^2)*(1-(yspan.^2/s^2))); % Sectional lift at each spanwise point (from equation of an ellipse)
LSec_zf=sqrt((LSecMax_zf^2)*(1-(yspan.^2/s^2)));

%LSec =  LSec_zf;
% Plotting lift distribution
figure
box on
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1.5);

plot(yspan,LSec,'-b','LineWidth',1.5)
hold on
plot(yspan,LSec_zf,'-r','LineWidth',1.5)
xlabel('Spanwise Position [m]')
ylabel('Sectional Lift [N/m]')
grid on

hold off

lift=trapz(yspan,LSec)
lift2=trapz(yspan,LSec_zf)

%INERTIA
%% Getting inertial load distribution of horizontal stabiliser - LINEAR

% Assume sectional weight varies proportionally to area i.e. chord - will
% give a trapezoidal weight distribution - area of this distribution must
% equal total weight

W=ULF*0.5*1230*g; % Weight of one horizontal stabiliser in N*****

semispan = bperp/2;
S1 = 6.96;
S2 = S1*(lambda)^2;
WingBox_Volume = 1/3*(S1+S2+(S1*S2)^(1/2))*semispan;

Area_Span = (S1*(lambda^2-2*lambda+1))*(yspan./semispan).^2 + (S1*(-2+2*lambda)).*(yspan./semispan) + S1;
WingW_Span = ((W)/WingBox_Volume).*Area_Span;
wSec = WingW_Span;
wSecMax=(2*W)/(s*(ctip/croot + 1))
wSecMin=(ctip/croot)*wSecMax
% % Plotting weight distribution
% figure
% box on
% set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold','LineWidth', 1.5);
% hold on
% grid on
% plot(yspan,WingW_Span,'r--',LineWidth=1.5)
% xlabel('Spanwise Position [m]')
% ylabel('Sectional Weight [N/m]')
% 
% hold off


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
%legend('Lift (MTOW)','Weight')
hold off


%% Getting shear force distribution of horizontal stabiliser

% Finding the length of segments dy
dy=yspan(2);

% Taking h.stab as a cantilever beam and taking cut to find shear force (in N/m):
for z=1:length(yspan) % z is spanwise coordinate but starting at the tip
    Shear(z) = dy*sum(LSec(length(LSec)-z+1:length(LSec))) - dy*sum(wSec(length(wSec)-z+1:length(wSec)));
    Shear_zf(z) = dy*sum(LSec_zf(length(LSec_zf)-z+1:length(LSec_zf))) - dy*sum(wSec(length(wSec)-z+1:length(wSec)));
end

% Flipping Shear to a coordinate system which starts from the root
Shear2=flip(Shear);
Shear2_zf=flip(Shear_zf);

% Plotting shear force distribution
figure
box on
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1.5);
hold on
plot(yspan,Shear2,'-b','LineWidth',1.5)
plot(yspan,Shear2_zf,'-r','LineWidth',1.5)

xlabel('Spanwise Position [m]')
ylabel('Shear Force [N]')
% legend('Load case 1 (MTOW)','Load Case 1 (MZFW)')
legend('MTOW Case','MZFW Case')
grid on
%ylim([0,10000]);
hold off


%% Getting bending moment distribution of horizontal stabiliser - LINEAR

% Cantilever beam analysis but convert lift and weight distributions into
% point forces at centroids of 1/4 ellipse and trapezium respectively

for i=1:length(yspan) % Spanwise coordinate starting from tip
    BeamLift = dy*sum(LSec(length(LSec)-i+1:length(LSec))); % Finding total lift up to section currently considering
    BeamLift_zf = dy*sum(LSec_zf(length(LSec_zf)-i+1:length(LSec_zf)));
    BeamWeight = dy*sum(wSec(length(wSec)-i+1:length(wSec))); % Finding total weight up to section currently considering
    currentw = wSec(length(wSec)-i+1); % Finding sectional weight at the cut (needed for finding centroid of trapezium so it finds its way into moment equation)
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
% BM2=abs(BM2);

BM2_zf=flip(BM_zf);
% BM2_zf=abs(BM2_zf);

% Plotting bending moment distribution
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

% T=La + nWb - M0, M0=0 for horizontal stabilisers since symmetrical airfoil

% Assumptions for location of aerodynamic centre, CG etc. - kept as
% variables so we can change later if we want/need to:
% AS A FRACTION OF CHORD LENGTH:
AC=0.25; % Aerodynamic centre
FA=0.45; % Flexural analysis
CG=0.391809984147464; % CG

% Finding torque distribution:
for i=1:length(yspan)
    SectionLiftMoment(i)=LSec(length(LSec)-i+1)*((FA-AC)*c(length(c)-i+1)); % Finding contribution of lift to torque of section
    SectionLiftMoment_zf(i)=LSec_zf(length(LSec)-i+1)*((FA-AC)*c(length(c)-i+1));
    SectionWeightMoment(i)=wSec(length(wSec)-i+1)*((CG-FA)*c(length(c)-i+1)); % Finding contribution of weight to torque of section
    TorqueSec(i)=SectionLiftMoment(i)+SectionWeightMoment(i); % Finding sectional torque
    Torque(i)= dy*sum(TorqueSec(length(TorqueSec)-i+1:length(TorqueSec))); % Finding torque at each section
    TorqueSec_zf(i)=SectionLiftMoment_zf(i)+SectionWeightMoment(i); % Finding sectional torque
    Torque_zf(i)= dy*sum(TorqueSec_zf(length(TorqueSec_zf)-i+1:length(TorqueSec_zf))); % Finding torque at each section
end

% Flipping to a coordinate system which starts from the root:
Torque2=flip(Torque);
Torque2_zf=flip(Torque_zf);

% Plotting torque distribution:
figure
box on
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 1.5);
hold on
grid on
plot(yspan,Torque2,'-b','LineWidth',1.5)
plot(yspan,Torque2_zf,'-r','LineWidth',1.5)
xlabel('Spanwise Position [m]')
ylabel('Torque [Nm]')
ylim([0,12000])
legend('MTOW Case','MZFW Case')
hold off



maxvalue = [max(abs([Shear,Shear2_zf])),max(abs([BM2,BM2_zf])),max(abs([Torque,Torque_zf]))]
