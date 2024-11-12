clc
clear

%Weights
W_o = 390057;        %MTOW
W_ng = 0.1*W_o;      %Nose gear weight carrying capability
W_mg = W_o - W_ng;     %Main gear weight carrying capability

MAC = 9.7; %in meters

Wo_lbs = W_o*2.20462; %Wo in lbs (errikos just being extremely annoying with values)
%This value is 85 9927 lbs

%Moments
x_wing = 37; % distance of wing leading edge 
x_cg = 37 +0.4*MAC;                             % distance of cg (estimate)
x_mg = 44.5;                               %x distance of main gears (estimate)
x_ng = (W_o*x_cg - W_mg*x_mg)/W_ng;     %x distance of nose gears



%Assumed values: tire pressure, number of u/c struts and wheels per strut

%Nose Gear: one struct, two wheels

%Main Gear: three struct, each with 4 wheels (twini tandem)

%Tire Sizing (statistical approach): main wheel
W_mw = W_mg/12; %Weight on each main wheel
D_mw = 1.63*W_mw^0.315; %Diameter of main wheel in INCH
Wth_mw =  0.1043*W_mw^0.48; %Width of main wheel in INCH

%Tire Sizing (statistical approach): Nose wheel
W_nw = W_ng/2; %Weight on each nose wheel
D_nw = 1.63*W_nw^0.315; %Diameter of nose wheel in INCH
Wth_nw = 0.1043*W_nw^0.48; %Width of nose wheel in INCH

%Tire Selection (From a table on Raymer's book), based on Max loading
%Three-part Name Tire: 47 x 18-18
R_r = 19.2; % Rolling Radius in inch

%Tire Pressure
A_p = 2.3*sqrt(Wth_mw*D_mw)*(D_mw/2 - R_r); %wheel contact area in inch^2
P_mw = W_mw/(A_p*0.00064516);

%W_o < 200 000 lbs so we need 2 wheels per main strut

%ACN - Aircraft Classification Number
%The follwoing values are generated from the COMFAA software

