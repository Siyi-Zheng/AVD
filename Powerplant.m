clear
clc

% Pitot Inlet


% Manufacturer Engine data
mdot = 2297; % takeoff mass flow rate in lbs/sec
D_fan = 104.7; % fan diamter in inches
L_engine = 169.7; % inches
fprintf("Fan diamter: %.3fm \n",D_fan*0.0254)
W = 6147.1;

%thrust required to meet T/W ratio
T_W = 0.29;

MTOW = 390000;

T_total = (T_W * MTOW) * 9.81;

T_req = (T_total / 4) ;



% Flight Condition
M_cruise = 0.83;

% unit conversion
D_fan = D_fan*0.0254; % inch to metre
L_engine = L_engine*0.0254; % inch to metre


%rubber engine corrections

%GEnx 1B75/P1 Engine Thrust
Thrust_to =29900; 

%propulsion intergration corrections

%inlet ram recovery factor
Cram = 1.35;

pr_ref = 1.0;

perc_thrustloss1 = Cram * (pr_ref - 0.98) * 100;

%Rubber engines scaling

Thrustloss = Thrust_to * (100 - perc_thrustloss1)/100;


% Capture Area 

% Empirical: according to raymer pg300 graph
A_c = mdot*0.025; % ft^2
D_c = sqrt(A_c/pi)*2*0.3048; % metre
fprintf("Capture area diameter (empirical): %.3fm \n",D_c)

% isentropic flow
M_engine = 0.4; % the speed in front of the engine fan
M_t = (M_cruise + M_engine)/2; % assumed to be half, according to Raymer and last year's report
A_ratio = (1/M_t*( (1+0.2*M_t^2)/1.2 )^3 ) / (1/M_engine*( (1+0.2*M_engine^2)/1.2 )^3);
D_ratio = sqrt(A_ratio);
D_th = D_ratio*D_fan; % metre
fprintf("Capture area diameter (isentropic): %.3fm \n",D_th)


% Diffuser Length (D_diff = 0.6-1 * D_fan)
diffuser_angle = 7; % in degree
L_diff = ((D_fan-D_th)/2) /tan(deg2rad(diffuser_angle));


% Highlight Diameter
D_1 = 0.92*D_fan; % factor 0.9-0.95


% Max macella diameter
D_max = 1.3*D_fan;


% Lip Geometry (a,b,c,d)
% (b+d) = L_diff*0.15 (factor from 0.15-0.2)
% b = 1.4*d (factor from 1.5-2)
d = (L_diff*0.15)/(1+1.4);
b = 1.4*d;
a = 2*b; % factor from 1.5-3
c = 4*d; % factor from 3-5
fprintf("a: %.3f, b: %.3f, c: %.3f, d: %.3f, \n",a,b,c,d)


% Nozzle
L_nozzle = 0.75*D_fan; % reference from book example
L_nacelle = L_diff + L_engine + L_nozzle;
L_nacelle/D_max % should be in range 2-3
D_nozzle = sqrt(0.75)*D_th; % factor 0.5-0.7 raymer, 0.5-0.75 cambridge book

% distance from fan to fan case 34.4cm


