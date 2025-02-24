function [M,N,S] = heavyframes(P , T , Q , R)

R = 6.34;

theta = (0:0.1:360) * pi/180;

% Tangential unit load case
M_T = (P*R/(2*pi)) * ((3*sin(theta)/2) + (pi - theta)*cos(theta) - 1);
N_T = (P/(2*pi)) * (sin(theta)/2 - (pi - theta)*cos(theta));
S_T = (P/(2*pi)) * ((pi - theta)*sin(theta) - (1 - cos(theta))/2);

% Radial unit load case
M_R = (Q*R/(2*pi)) * (cos(theta)/2 + (pi - theta)*sin(theta) + 1);
N_R = (Q/(2*pi)) * (3*cos(theta)/2 + (pi - theta)*sin(theta));
S_R = (Q/(2*pi)) * ((pi - theta)*cos(theta) - sin(theta)/2);

% Moment unit load case
M_M = (T/(2*pi)) * (pi - 2*sin(theta) - theta);
N_M = (T/(2*pi*R)) * (3*cos(theta)/2 + (pi - theta)*sin(theta));
S_M = (T/(2*pi*R)) * (1 + 2*cos(theta));

%superpose final results
M = M_T + M_R + M_M;
N = N_T + N_R + N_M;
S = S_T + S_R + S_M;


end