clear all
clc

W_landing = 0.85 * 353385 * 9.81; % landing weight (N)
beta = 11; % tipback angle (deg)
V_v = 3.05; % vertical speed at landing (m/s)
n = 3; % landing load factor
F_gear = W_landing * (n-1); % vertical main gear force
D_gear = 0.8 * W_landing; % horizontal main gear force (drag)
% horizontal and vertical forces transmitted to airframe
H = F_gear * sind(beta) - D_gear * cosd(beta);
V = F_gear * cosd(beta) + D_gear * sind(beta);