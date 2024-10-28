clear
clc

% Pitot Inlet


% Manufacturer Engine data
mdot = 2559; % takeoff mass flow rate in lbs/sec
D_fan = 111.1; % fan diamter in inches
fprintf("Fan diamter: %.3fm \n",D_fan*0.0254)
M = 0.83;

%propulsion intergration corrections

%inlet ram recovery factor
Cram = 1.35;

pr_ref = 1.0;

perc_thrustloss = Cram * (pr_ref - 0.98) * 100;


% Capture Area 

% according to raymer pg300 graph
A_c = mdot*0.025; % ft^2
D_c = sqrt(A_c/pi)*2*0.3048; % metre
fprintf("Capture area diameter (empirical): %.3fm \n",D_c)

% isentropic flow
M_engine = 0.4; % the speed in front of the engine fan
M_cruise = 0.83;
M_t = (M_cruise + M_engine)/2; % assumed to be half, according to Raymer and last year's report

A_ratio = (1/M_t*( (1+0.2*M_t^2)/1.2 )^3 ) / (1/M_engine*( (1+0.2*M_engine^2)/1.2 )^3);
D_ratio = sqrt(A_ratio);
D_t = D_ratio*D_fan*0.0254; % metre
fprintf("Capture area diameter (isentropic): %.3fm \n",D_t)

% Nozzle

