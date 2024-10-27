clear
clc 
close all

%wings optimisation

W_S_takeoff= 7637;
Sref= 553.6;
MTOW= 430411;
fuel_weight= 198928;
T_W_TO= 0.2898;

MWF=[0.970
0.985
0.611
0.990
0.985
0.988
0.984
0.990 
0.995
]';

weight_average_cruise = (MWF(1)*MWF(2)*0.8055)*MTOW*9.81;
q_cruise= 0.5*0.3016*(0.83 * 295)^2;

design_cl_average= (1/q_cruise) * (weight_average_cruise/Sref);



