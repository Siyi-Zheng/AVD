clear all
clc

W_landing = 0.85 * 353385 * 9.81; % landing weight (N)
V_v = 3.05; % vertical speed at landing (m/s)
n = 3; % landing load factor
F_gear = (W_landing * n) / 4; % divide between the 4 landing gear