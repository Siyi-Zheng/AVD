clc
clear
%Weights
W_o = 390057;        %MTOW
W_ng = 0.12*W_o;      %Nose gear weight carrying capability
W_mg = W_o - W_ng;     %Main gear weight carrying capability

Wo_lbs = W_o*2.20462; %Wo in lbs (errikos just being extremely annoying with values)
%This value is 85 9927 lbs

%Moments
x_cg = 37;                              %x distance of cg (estimate)
x_ng = 5;                               %x distance of nose gears (estimate)
x_mg = (W_o*x_cg - W_mg*x_ng)/W_mg;     %x distance of main gears

%Assumed values: tire pressure, number of u/c struts and wheels per strut
%W_o < 200 000 lbs so we need 2 wheels per main strut

%ACN - Aircraft Classification Number
%The follwoing values are generated from the COMFAA software

