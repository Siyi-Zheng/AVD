clear all
clc
close all

%Section - Structural Idealisation

%Wing Box Plot
%Spars at 0.15 & 0.75.
fspar_pos = 0.12;
rspar_pos = 0.7;
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

upper =  [0.000000  0.000000
  0.002000  0.010770
  0.005000  0.016580
  0.010000  0.022400
  0.020000  0.029600
  0.030000  0.034600
  0.040000  0.038300
  0.050000  0.041400
  0.060000  0.044000
  0.070000  0.046300
  0.080000  0.048400
  0.090000  0.050200
  0.100000  0.051900
  0.110000  0.053500
  0.120000  0.054900
  0.130000  0.056200
  0.140000  0.057400
  0.150000  0.058600
  0.160000  0.059700
  0.170000  0.060700
  0.180000  0.061600
  0.190000  0.062500
  0.200000  0.063300
  0.210000  0.064100
  0.220000  0.064800
  0.230000  0.065400
  0.240000  0.066000
  0.250000  0.066500
  0.260000  0.067000
  0.270000  0.067500
  0.280000  0.067900
  0.290000  0.068300
  0.300000  0.068600
  0.310000  0.068900
  0.320000  0.069200
  0.330000  0.069400
  0.340000  0.069600
  0.350000  0.069700
  0.360000  0.069800
  0.370000  0.069900
  0.380000  0.069900
  0.390000  0.069900
  0.400000  0.069900
  0.410000  0.069800
  0.420000  0.069700
  0.430000  0.069600
  0.440000  0.069500
  0.450000  0.069300
  0.460000  0.069100
  0.470000  0.068900
  0.480000  0.068600
  0.490000  0.068300
  0.500000  0.068000
  0.510000  0.067600
  0.520000  0.067200
  0.530000  0.066800
  0.540000  0.066300
  0.550000  0.065800
  0.560000  0.065300
  0.570000  0.064700
  0.580000  0.064100
  0.590000  0.063500
  0.600000  0.062800
  0.610000  0.062100
  0.620000  0.061300
  0.630000  0.060500
  0.640000  0.059700
  0.650000  0.058800
  0.660000  0.057900
  0.670000  0.056900
  0.680000  0.055900
  0.690000  0.054800
  0.700000  0.053700
  0.710000  0.052500
  0.720000  0.051300
  0.730000  0.050000
  0.740000  0.048700
  0.750000  0.047300
  0.760000  0.045800
  0.770000  0.044300
  0.780000  0.042700
  0.790000  0.041100
  0.800000  0.039400
  0.810000  0.037600
  0.820000  0.035800
  0.830000  0.033900
  0.840000  0.031900
  0.850000  0.029900
  0.860000  0.027800
  0.870000  0.025600
  0.880000  0.023400
  0.890000  0.021100
  0.900000  0.018700
  0.910000  0.016200
  0.920000  0.013700
  0.930000  0.011100
  0.940000  0.008400
  0.950000  0.005600
  0.960000  0.002700
  0.970000 -0.000200
  0.980000 -0.003200
  0.990000 -0.006300
  1.000000 -0.009500];

lower = [0.000000  0.000000
  0.002000 -0.010770
  0.005000 -0.016580
  0.010000 -0.022400
  0.020000 -0.029600
  0.030000 -0.034500
  0.040000 -0.038200
  0.050000 -0.041300
  0.060000 -0.043900
  0.070000 -0.046200
  0.080000 -0.048300
  0.090000 -0.050100
  0.100000 -0.051800
  0.110000 -0.053400
  0.120000 -0.054900
  0.130000 -0.056200
  0.140000 -0.057400
  0.150000 -0.058600
  0.160000 -0.059700
  0.170000 -0.060700
  0.180000 -0.061600
  0.190000 -0.062500
  0.200000 -0.063300
  0.210000 -0.064100
  0.220000 -0.064800
  0.230000 -0.065500
  0.240000 -0.066100
  0.250000 -0.066700
  0.260000 -0.067200
  0.270000 -0.067700
  0.280000 -0.068100
  0.290000 -0.068500
  0.300000 -0.068800
  0.310000 -0.069100
  0.320000 -0.069300
  0.330000 -0.069500
  0.340000 -0.069600
  0.350000 -0.069700
  0.360000 -0.069700
  0.370000 -0.069700
  0.380000 -0.069600
  0.390000 -0.069500
  0.400000 -0.069300
  0.410000 -0.069100
  0.420000 -0.068800
  0.430000 -0.068500
  0.440000 -0.068100
  0.450000 -0.067700
  0.460000 -0.067200
  0.470000 -0.066700
  0.480000 -0.066100
  0.490000 -0.065400
  0.500000 -0.064600
  0.510000 -0.063700
  0.520000 -0.062700
  0.530000 -0.061600
  0.540000 -0.060400
  0.550000 -0.059100
  0.560000 -0.057700
  0.570000 -0.056200
  0.580000 -0.054600
  0.590000 -0.052900
  0.600000 -0.051100
  0.610000 -0.049200
  0.620000 -0.047300
  0.630000 -0.045300
  0.640000 -0.043300
  0.650000 -0.041200
  0.660000 -0.039100
  0.670000 -0.037000
  0.680000 -0.034800
  0.690000 -0.032600
  0.700000 -0.030400
  0.710000 -0.028200
  0.720000 -0.026000
  0.730000 -0.023800
  0.740000 -0.021600
  0.750000 -0.019400
  0.760000 -0.017300
  0.770000 -0.015200
  0.780000 -0.013200
  0.790000 -0.011300
  0.800000 -0.009500
  0.810000 -0.007900
  0.820000 -0.006400
  0.830000 -0.005000
  0.840000 -0.003800
  0.850000 -0.002800
  0.860000 -0.002000
  0.870000 -0.001400
  0.880000 -0.001000
  0.890000 -0.000800
  0.900000 -0.000900
  0.910000 -0.001200
  0.920000 -0.001700
  0.930000 -0.002500
  0.940000 -0.003600
  0.950000 -0.005000
  0.960000 -0.006700
  0.970000 -0.008700
  0.980000 -0.011000
  0.990000 -0.013600
  1.000000 -0.016500];

height = upper - lower; % get height of spar along x/c
centroid = upper + lower; % get centroid of spar along x/c
h1 = upper(:,1);
h2 = height(:,2);
c = centroid(:,2);
final = cumtrapz(h1, h2); % integrate the height to find the mean
final2 = cumtrapz(h1, c); % integrate the centroid to find the mean
spar_height = (h2(73) + h2(15)) / 2; % mean value theorem
box_centroid = (final2(73) - final2(15)) / (0.7 - 0.12); % mean value theorem
box_lower = box_centroid - spar_height/2; % find the lower end of the wing box

plot(NACA(:,1),NACA(:,2),'k',LineWidth=1.1)
axis equal 
grid on
xlabel('Chordwise position (x/c)')
ylabel('Streamwise position (y/c)')
xlim([0,max(NACA(:,1))])
ylim([-0.1 0.1])
hold on

rectangle('Position',[fspar_pos,box_lower,rspar_pos-fspar_pos,spar_height],'EdgeColor','r','LineWidth',1.3)

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
figure(1)
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

L0 = (4*3.75*9.81*353000)/(pi*semispan*2);

for i = 1:length(yspan)
    if yspan(i) >= 0 %fuselage_diam/2
        Spanwise_Lift(i) = L0*(1-(yspan(i)/(semispan))^2)^(1/2);
    else
        Spanwise_Lift(i) = 0;
    end
end

%Weight of ONE Wing in N
Wing_Weight = 9.81*27100/2;
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

% Flight speeds for load cases
VA = 250; % Maneuver speed (m/s)
VD = 320; % Dive speed (m/s)

% Symmetric Flight Lift Distribution
L0_VA = (4 * 3.75 * 9.81 * 353000) / (pi * semispan * 2); % Ultimate load factor at VA
L0_VD = (4 * 3.0 * 9.81 * 353000) / (pi * semispan * 2); % Reduced load factor at VD

Spanwise_Lift_VA = zeros(size(yspan));
Spanwise_Lift_VD = zeros(size(yspan));
for i = 1:length(yspan)
    Spanwise_Lift_VA(i) = L0_VA * sqrt(1 - (yspan(i) / semispan)^2);
    Spanwise_Lift_VD(i) = L0_VD * sqrt(1 - (yspan(i) / semispan)^2);
end

% One Engine Inoperative (OEI) Load
engine_thrust = 150000; % Thrust of one engine in Newtons
engine_offset = 10; % Distance from CG to engine in meters
vertical_stab_offset = 30; % Distance from CG to vertical stabilizer aerodynamic center in meters

yaw_moment = engine_thrust * engine_offset; % Yawing moment due to OEI
stab_force = yaw_moment / vertical_stab_offset; % Horizontal stabilizer force

% Apply stabilizer force as a horizontal load
OEI_Load = zeros(size(yspan));
stab_position = semispan + 5; % Approximate spanwise location of the tailplane
for i = 1:length(yspan)
    if abs(yspan(i) - stab_position) <= 1 % Apply load over a 2-meter span
        OEI_Load(i) = stab_force / 2; % Uniform distribution over 2 meters
    end
end


%ADD IN LANDING GEAR WEIGHT
for i = 1:length(yspan)
    if yspan(i) < 5.6
        WingW_Span(i) = WingW_Span(i) + Wing_LandingG*3.75/(5.6);
        Gear_Load_Span(i) = Gear_Load_Span(i) + Landing_GearLoad/5.6;
    end
end

% ADD ENGINE WEIGHT TO THE SPANWISE DISTRIBUTION
engine_mass = 6550; % Mass of one engine in kilograms
engine_load = engine_mass * 9.81; % Weight of one engine in Newtons
engine_positions = [9.75, 21.1]; % Spanwise positions of the engines in meters

for i = 1:length(yspan)
    % Check if within the 1-meter influence region (0.5 m on each side) of the first engine
    if yspan(i) >= engine_positions(1) - 0.5 && yspan(i) <= engine_positions(1) + 0.5
        WingW_Span(i) = WingW_Span(i) + engine_load * 3.75 / 1.0; % Distribute over 1 meter
        Gear_Load_Span(i) = Gear_Load_Span(i) + engine_load / 1.0; % Distribute over 1 meter
    end
    % Check if within the 1-meter influence region (0.5 m on each side) of the second engine
    if yspan(i) >= engine_positions(2) - 0.5 && yspan(i) <= engine_positions(2) + 0.5
        WingW_Span(i) = WingW_Span(i) + engine_load * 3.75 / 1.0; % Distribute over 1 meter
        Gear_Load_Span(i) = Gear_Load_Span(i) + engine_load / 1.0; % Distribute over 1 meter
    end
end



Landing_Spanwise_Lift = Spanwise_Lift.*0.85;
LandingW_Span = -WingW_Span-(0.85.*FuelW_Span-0.15*WingW_Span);
figure(2)
plot(yspan-1.305,Spanwise_Lift,'b')
hold on
plot(yspan, Spanwise_Lift_VA, 'b', 'LineWidth', 1.5);
plot(yspan, Spanwise_Lift_VD, 'r', 'LineWidth', 1.5);
plot(yspan, OEI_Load, 'g', 'LineWidth', 1.5);
plot(yspan-1.305,Landing_Spanwise_Lift,'b:')
plot(yspan-1.305,Spanwise_Lift.*(24368/48725),'b--')
plot(yspan-1.305,-WingW_Span-FuelW_Span,'r')
plot(yspan-1.305,LandingW_Span,'r:')
plot(yspan-1.305,-WingW_Span,'r--')
plot(yspan-1.305,Gear_Load_Span,'g:',LineWidth=1)
legend('Lift Force MTOW','Symmetric Flight - VA', 'Symmetric Flight - VD', 'OEI Load''Lift Force MLW','Lift Force EFW','Wing Sectional Weight MTOW','Wing Sectional Weight MLW','Wing Sectional Weight EFW','Landing Gear Load')
xlabel('Spanwise Position (m)')
ylabel('Sectional Load (N/m)')
grid 
xlim([0,33])
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
        MTOW_WShear(n-i) = MTOW_WShear(n-i) + Wing_LandingG*3.75;
        EFW_WShear(n-i) = EFW_WShear(n-i) + Wing_LandingG*3.75;
        MLW_WShear(n-1) = MLW_WShear(n-i) + Wing_LandingG*3.75;
    end
end
figure
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
xlim([0,33])
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
figure
plot(yspan-1.305,MTOW_BM,'b')
hold on
plot(yspan-1.305,EFW_BM,'r')
plot(yspan-1.305,MLW_BM,'g')
xlabel('Span Position (m)')
ylabel('Bending Moment (Nm)')
legend('LC 1: MTOW','LC 1: EFW','LC 3: MLW')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth',1.5);
grid on
xlim([0,33])
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
figure
plot(yspan-1.305,Torque_MTOW,'b')
hold on
plot(yspan-1.305,Torque_EW,'r')
plot(yspan-1.305,Torque_landing,'g')
grid on
legend('LC 1: MTOW','LC 1: EFW','LC 3: MLW')
xlabel('Spanwise Position (m)')
ylabel('Torque (Nm)')
xlim([0,33])
set(findobj(gcf, 'type', 'axes'),'FontSize', 13, 'FontWeight', 'Bold', 'LineWidth', 1);
set(findobj(gcf, 'type', 'line'), 'LineWidth', 1.5);
xlabel(get(get(gca,'XLabel'),'String'),'Interpreter','latex');
ylabel(get(get(gca,'YLabel'),'String'),'Interpreter','latex');
lg = legend;
set(lg, 'Interpreter','latex');
hold off
