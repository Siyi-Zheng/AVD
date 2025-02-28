function [M,N,S,theta] = heavyframes(P , T , Q )

R = 3.17;

theta = (0:0.1:360) * pi/180;

% Tangential unit load case
M_P = ((P*R)/(2*pi)) .* ((3*sin(theta)./2) + (pi - theta).*(cos(theta) - 1));
N_P = (P/(2*pi)) .* (sin(theta)./2 - (pi - theta).*cos(theta));
S_P = (P/(2*pi)) .* ((pi - theta).*sin(theta) - 1 - (cos(theta)./2));

% Radial unit load case
M_Q = ((Q*R)/(2*pi)) .* (cos(theta)./2 - (pi - theta).*sin(theta) + 1);
N_Q = (Q/(2*pi)) .* (3.*cos(theta)./2 + (pi - theta).*sin(theta));
S_Q = (Q/(2*pi)) .* ((pi - theta).*cos(theta) - sin(theta)./2);

% Moment unit load case
M_T = (T/(2*pi)) .* (pi - 2.*sin(theta) - theta);
N_T = (T/(2*pi*R)) .* (3.*cos(theta)/2 + (pi - theta).*sin(theta));
S_T = (T/(2*pi*R)) .* (1 + 2.*cos(theta));

%superpose final results
M = M_T + M_Q + M_P;
N = N_T + N_Q + N_P;
S = S_T + S_Q + S_P;


end