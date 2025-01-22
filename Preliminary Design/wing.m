clear all
clc
close all

%Section - Structural Idealisation

%Wing Box Plot
%Spars at 0.15 & 0.75.
fspar_pos = 0.15;
rspar_pos = 0.7;
spar_height = (0.074221+0.062372)/2;
fuselage_diam = 6.34;

NACA = [1.000  -.0104
0.99    -.0071
0.98    -.0039
0.97    -.0009
0.95    .0049
0.92    .0131
0.9      .0181
0.87    .0251
0.85    .0294
0.82    .0353
0.800  .0389
0.77    .0439
0.75    .0469
0.72    .0509
0.7      .0533
0.68    .0555
0.65    .0585
0.62    .0610
0.6      .0625
0.57    .0645
0.55    .0656
0.53    .0666
0.50    .0678
0.48    .0684
0.45    .0692
0.43    .0695
0.4      .0697
0.38    .0698
0.35    .0696
0.33    .0692
0.3      .0685
0.27    .0673
0.25    .0664
0.22    .0646
0.2      .0632
0.17    .0606
0.15    .0585
0.12    .0548
0.1      .0518
0.07    .0462
0.05    .0411
0.04    .0381
0.03    .0343
0.02    .0293
0.01    .0219
0.005  .0158
0.002  .0095
0.00    .00
0.002  -.0093
.005    -.016
0.01    -.0221
0.02    -.0295
0.03    -.0344
0.04    -.0381
0.05    -.0412
0.07    -.0462
0.1      -.0517
0.12    -.0547
0.15    -.0585
0.17    -.0606
0.20    -.0633
0.22    -.0647
0.25    -.0666
0.28    -.068
0.3      -.0687
0.32    -.0692
0.35    -.0696
0.37    -.0696
0.4      -.0692
0.42    -.0688
0.45    -.0676
0.48    -.0657
0.5      -.0644
0.53    -.0614
0.55    -.0588
0.58    -.0543
0.6      -.0509
0.63    -.0451
0.65    -.041
0.68    -.0346
0.70    -.0302
0.73    -.0235
0.75    -.0192
0.77    -.0150
0.80    -.0093
0.83    -.0048
0.85    -.0024
0.87    -.0013
0.89   -.0008
0.92    -.0016
0.94    -.0035
0.95    -.0049
0.96    -.0066
0.97    -.0085
0.98    -.0109
0.99    -.0137
1.0      -.0163];

plot(NACA(:,1),NACA(:,2),'k',LineWidth=1.1)
axis equal 
grid on
xlabel('X-pos x/c')
ylabel('Y-pos y/c')
xlim([0,max(NACA(:,1))])
hold on

rectangle('Position',[fspar_pos,-0.025555+0.006,rspar_pos-fspar_pos,spar_height],'EdgeColor','r','LineWidth',1.3)

hold off

%Wing Plot
xpos_wing = 0.0;
sweep = 30;
AR = 8.77;
Sref = 482; %both wings
lambda = 0.25; %taper ratio 

semispan = 0.5*(AR*Sref)^(1/2);
wing_CG = [13.03,1.3000];

n = 500;                            %No. of Stations
yspan = linspace(0,semispan,n);

wing_root = Sref/(semispan*(1+lambda));

root_sparheight = spar_height.*wing_root;
spar_height_span = (lambda-1)*root_sparheight.*(yspan./semispan)+root_sparheight;

yleading = linspace(0,semispan,100);
xleading = tan(sweep*pi/180).*yleading + xpos_wing;

xroot = linspace(0,wing_root,100) + xpos_wing;
yroot = ones([100,1]).*0;

xtip = linspace(0,wing_root*lambda,100) + max(xleading);
ytip = ones([100,1]).*max(yleading);

xtrailing = linspace(xroot(length(xroot)),xtip(length(xtip)),100);
ytrailing = ((ytip(length(ytip))-yroot(length(yroot))))/(xtip(length(xtip))-xroot(length(xroot))).*(xtrailing) + (yroot(length(yroot)) - xroot(length(xroot))*((ytip(length(ytip))-yroot(length(yroot))))/(xtip(length(xtip))-xroot(length(xroot))));

flexural_axis_x = linspace(wing_root*(rspar_pos+fspar_pos)/2,min(xtip)+lambda*wing_root*(rspar_pos+fspar_pos)/2,n);
skin_COM_x = linspace(0.4167*wing_root,min(xtip)+0.4167*lambda*wing_root,n);
spar_COM_x = linspace(wing_root*(rspar_pos+fspar_pos)/2,min(xtip)+lambda*wing_root*(rspar_pos+fspar_pos)/2,n);
centreline_x = (xleading+xtrailing)./2;
fuel_COM_x = linspace(wing_root*(rspar_pos+fspar_pos)/2,min(xtip)+lambda*wing_root*(rspar_pos+fspar_pos)/2,n);
rspar_x = linspace(wing_root*rspar_pos,min(xtip)+wing_root*lambda*rspar_pos,n);

qchord_x = (xleading+centreline_x)./2;
centreline_x = interp(centreline_x,n/100);
qchord_x = interp(qchord_x,n/100);

%Wing Plots
% - Right Wing
plot(xleading,yleading,'k',LineWidth=1)
hold on
plot(xtip,ytip,'k',LineWidth=1)
plot(xtrailing,ytrailing,'k',LineWidth=1)
line([fspar_pos*wing_root+xpos_wing,min(xtip)+fspar_pos*(max(xtip)-min(xtip))],[0,semispan],'Color','red','LineWidth',1.3)
line([rspar_pos*wing_root+xpos_wing,min(xtip)+rspar_pos*(max(xtip)-min(xtip))],[0,semispan],'Color','red','LineWidth',1.3)
plot(4.89,1.8,'bx',LineWidth=1.5)
plot(flexural_axis_x,yspan,'b--')
%plot(centreline_x,yspan,'g--')
%plot(qchord_x,yspan,'g--')
%lot(skin_COM_x,yspan,'g--')
%plot(wing_CG(1),wing_CG(2),'ko')#
grid on
axis equal
hold off

%Secton Wing Loads
%L0 uses Empty Fuel Weight as the Worst Case:

L0 = (2*2.5*9.81*162976)/(pi*semispan*2);

for i = 1:length(yspan)
    if yspan(i) >= 0 %fuselage_diam/2
        Spanwise_Lift(i) = L0*(1-(yspan(i)/(semispan))^2)^(1/2);
    else
        Spanwise_Lift(i) = 0;
    end
end

%Weight of ONE Wing in N
Wing_Weight = 9.81*4968/2;
%Weight of ONE Wing fuel in N
Wing_MaxFuel = 9.81*(0.413*(162976-27100))/2;
%Weight of ONE landing gear in N
Wing_LandingG = 9.81*(16650)/2;
%Load Case 3 Gear Load in N
Landing_GearLoad = 195030/2;

S1 = wing_root*(rspar_pos-fspar_pos)*spar_height*wing_root;
S2 = S1*(lambda)^2;
WingBox_Volume = 1/3*(S1+S2+(S1*S2)^(1/2))*semispan;

Area_Span = (S1*(lambda^2-2*lambda+1)).*(yspan./semispan).^2 + (S1*(-2+2*lambda)).*(yspan./semispan) + S1;
WingW_Span = ((Wing_Weight*3.75)/WingBox_Volume).*Area_Span;
Gear_Load_Span = zeros([1,length(yspan)]);
FuelW_Span = zeros([1,length(yspan)]);

for i = 1:length(yspan)
    if yspan(i) <= 0.9*semispan
        FuelW_Span(i) = ((Wing_MaxFuel*2.5)/(0.9*WingBox_Volume))*Area_Span(i);
    end
end

%ADD IN LANDING GEAR WEIGHT
for i = 1:length(yspan)
    if yspan(i) < 1.8
        WingW_Span(i) = WingW_Span(i) + Wing_LandingG*2.5/(1.8);
        Gear_Load_Span(i) = Gear_Load_Span(i) + Landing_GearLoad/1.8;
    end
end


Landing_Spanwise_Lift = Spanwise_Lift.*0.85;
LandingW_Span = -WingW_Span-(0.85.*FuelW_Span-0.15*WingW_Span);
plot(yspan-1.305,Spanwise_Lift,'b')
hold on
plot(yspan-1.305,Landing_Spanwise_Lift,'b:')
plot(yspan-1.305,Spanwise_Lift.*(24368/48725),'b--')
plot(yspan-1.305,-WingW_Span-FuelW_Span,'r')
plot(yspan-1.305,LandingW_Span,'r:')
plot(yspan-1.305,-WingW_Span,'r--')
plot(yspan-1.305,Gear_Load_Span,'g:',LineWidth=1)
legend('Lift Force MTOW','Lift Force MLW','Lift Force EFW','Wing Sectional Weight MTOW','Wing Sectional Weight MLW','Wing Sectional Weight EFW','Landing Gear Load')
xlabel('Spanwise Position (m)')
ylabel('Sectional Load (N/m)')
grid 
xlim([0,13])
set(findobj(gcf, 'type', 'axes'),'FontSize', 13, 'FontWeight', 'Bold', 'LineWidth', 1);
set(findobj(gcf, 'type', 'line'), 'LineWidth', 1.5);
xlabel(get(get(gca,'XLabel'),'String'),'Interpreter','latex');
ylabel(get(get(gca,'YLabel'),'String'),'Interpreter','latex');
lg = legend;
set(lg, 'Interpreter','latex');
hold off


%Derivation of Shear Force Diagrams

MTOW_LShear = cumtrapz(flip(yspan),flip(Spanwise_Lift));
EFW_LShear = cumtrapz(flip(yspan),flip(Spanwise_Lift.*(27100/162976)));

MTOW_WShear = cumtrapz(flip(yspan),flip(-WingW_Span-FuelW_Span));
EFW_WShear = cumtrapz(flip(yspan),flip(-WingW_Span));

MLW_LShear = cumtrapz(flip(yspan),flip(Landing_Spanwise_Lift));
MLW_WShear = cumtrapz(flip(yspan),flip(LandingW_Span));
MLW_GearShear = cumtrapz(flip(yspan),flip(Gear_Load_Span));

for i = 1:n
    if yspan(i) <= 1.8
        MTOW_WShear(n-i) = MTOW_WShear(n-i) + Wing_LandingG*2.5;
        EFW_WShear(n-i) = EFW_WShear(n-i) + Wing_LandingG*2.5;
        MLW_WShear(n-1) = MLW_WShear(n-i) + Wing_LandingG*2.5;
    end
end
plot(flip(yspan),-MTOW_LShear,'b')
hold on
plot(flip(yspan),-MLW_LShear,'b:')
plot(flip(yspan),-EFW_LShear,'b--')
plot(flip(yspan),-MTOW_WShear,'r')
plot(flip(yspan),-MLW_WShear,'r:')
plot(flip(yspan),-EFW_WShear,'r--')
plot(flip(yspan),-MLW_GearShear,'g:')
legend('Lift Shear MTOW (N)','Lift Shear MLW (N)','Lift Shear EFW (N)','Wing Weight Shear MTOW (N)','Wing Weight Shear MLW (N)','Wing Weight Shear EFW (N)','Gear Load Shear',Location='northeast')
title('Load Case 1 MTOW vs. EFW Shear Diagrams')
hold off
plot(flip(yspan)-1.305,-(MTOW_LShear+MTOW_WShear),'b')
hold on
plot(flip(yspan)-1.305,-(EFW_LShear+EFW_WShear),'r')
plot(flip(yspan)-1.305,-(MLW_LShear+MLW_WShear+MLW_GearShear),'g')
grid on
legend('LC 1: MTOW','LC 1: EFW','LC 3: EFW',Location='northeast')
ylabel('Shear Force (N)')
xlabel('Spanwise Position (m)')
xlim([0,13])
set(findobj(gcf, 'type', 'axes'),'FontSize', 13, 'FontWeight', 'Bold', 'LineWidth', 1);
set(findobj(gcf, 'type', 'line'), 'LineWidth', 1.5);
xlabel(get(get(gca,'XLabel'),'String'),'Interpreter','latex');
ylabel(get(get(gca,'YLabel'),'String'),'Interpreter','latex');
lg = legend;
set(lg, 'Interpreter','latex');
hold off

%Bending Moments
MTOW_BM = flip(cumtrapz((yspan),-(MTOW_LShear+MTOW_WShear)));
EFW_BM = flip(cumtrapz((yspan),-(EFW_LShear+EFW_WShear)));
MLW_BM = flip(cumtrapz((yspan),-(MLW_LShear+MLW_WShear+MLW_GearShear)));
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth',1.5);
plot(yspan-1.305,MTOW_BM,'b')
hold on
plot(yspan-1.305,EFW_BM,'r')
plot(yspan-1.305,MLW_BM,'g')
xlabel('Span Position (m)')
ylabel('Bending Moment (Nm)')
legend('LC 1: MTOW','LC 1: EFW','LC 3: MLW')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth',1.5);
grid on
xlim([0,13])
set(findobj(gcf, 'type', 'axes'),'FontSize', 13, 'FontWeight', 'Bold', 'LineWidth', 1);
set(findobj(gcf, 'type', 'line'), 'LineWidth', 1.5);
xlabel(get(get(gca,'XLabel'),'String'),'Interpreter','latex');
ylabel(get(get(gca,'YLabel'),'String'),'Interpreter','latex');
lg = legend;
set(lg, 'Interpreter','latex');
hold off

%Torque 
chord_x = (wing_root*(lambda-1))/(semispan).*yspan + wing_root;
CM_EW = spar_COM_x.*0.8 + skin_COM_x.*0.2;
CM_MTOW = (CM_EW.*Wing_Weight + fuel_COM_x.*(Wing_MaxFuel-Wing_Weight))/(Wing_MaxFuel);


M0_cruise = (1/2*(0.287)*(0.81*295)^2).*(chord_x.^2).*(-0.081);
M0_landing = (1/2*(1.23)*(91.8*1.15)^2).*(chord_x.^2).*(-0.081);

%Torque_EW = (Spanwise_Lift.*(24368/48725)).*(-qchord_x+flexural_axis_x) + WingW_Span.*(-CM_EW+flexural_axis_x)-M0_cruise;
%Torque_MTOW = ().*((Spanwise_Lift).*(-qchord_x+flexural_axis_x) + (WingW_Span+FuelW_Span).*(-CM_MTOW+flexural_axis_x))-M0_cruise;
%Torque_landing= (Spanwise_Lift.*0.85).*(-qchord_x+flexural_axis_x) + (LandingW_Span).*(-CM_MTOW+flexural_axis_x)-M0_cruise;

for i = 1:n-1
    Torque_MTOW_calc(i) = (Spanwise_Lift(i))*(-qchord_x(i)+flexural_axis_x(i))*(yspan(i+1)-yspan(i))+(WingW_Span(i)+FuelW_Span(i))*(-CM_MTOW(i)+flexural_axis_x(i))*(yspan(i+1)-yspan(i))-M0_cruise(i)*(yspan(i+1)-yspan(i));
    Torque_EW_calc(i) = (Spanwise_Lift(i)*(24368/48725)).*(-qchord_x(i)+flexural_axis_x(i))*(yspan(i+1)-yspan(i)) + WingW_Span(i)*(-CM_EW(i)+flexural_axis_x(i))*(yspan(i+1)-yspan(i))-M0_cruise(i)*(yspan(i+1)-yspan(i));
    Torque_landing_calc(i) =  (Spanwise_Lift(i)*0.85)*(-qchord_x(i)+flexural_axis_x(i))*(yspan(i+1)-yspan(i))+(LandingW_Span(i))*(-CM_MTOW(i)+flexural_axis_x(i))*(yspan(i+1)-yspan(i))-(rspar_x(i)-flexural_axis_x(i))*(Gear_Load_Span(i))*(yspan(i+1)-yspan(i))-M0_landing(i)*(yspan(i+1)-yspan(i));
end
for i = 1:n-1
    temp = Torque_MTOW_calc(i:n-1); 
    temp2 = Torque_EW_calc(i:n-1);
    temp3 = Torque_landing_calc(i:n-1);
    Torque_MTOW(i) = sum(temp);
    Torque_EW(i) = sum(temp2);
    Torque_landing(i) = sum(temp3);
end

Torque_MTOW(n) = 0;
Torque_EW(n) = 0;
Torque_landing(n) = 0;
% Torque_MTOW = Torque_MTOW - M0_cruise;
% Torque_EW = Torque_EW - M0_cruise;
% Torque_landing = Torque_landing - M0_landing;
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth',1.5);
plot(yspan-1.305,Torque_MTOW,'b')
hold on
plot(yspan-1.305,Torque_EW,'r')
plot(yspan-1.305,Torque_landing,'g')
grid on
legend('LC 1: MTOW','LC 1: EFW','LC 3: MLW')
xlabel('Spanwise Position (m)')
ylabel('Torque (Nm)')
xlim([0,13])
set(findobj(gcf, 'type', 'axes'),'FontSize', 13, 'FontWeight', 'Bold', 'LineWidth', 1);
set(findobj(gcf, 'type', 'line'), 'LineWidth', 1.5);
xlabel(get(get(gca,'XLabel'),'String'),'Interpreter','latex');
ylabel(get(get(gca,'YLabel'),'String'),'Interpreter','latex');
lg = legend;
set(lg, 'Interpreter','latex');
hold off
