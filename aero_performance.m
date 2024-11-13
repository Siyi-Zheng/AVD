clear
clc

AR = 8.77;
Cla = 7.16;
d = 6.34;
b = 65;
Sexp = 440;
Sref = 482;
sweep = 24; % at max thickness pt.
le_sweep = 30;

CLa = [];
M_list1 = [];
M_list2 = [];
for M = 0:0.01:2

    beta = (1 - M ^ 2) ^ 0.5;
    eta = beta * Cla / (2 * pi);
    F = 1.07 * (1 + d / b) ^ 2;

    CLa_subsonic = 2 * pi * AR * Sexp * F /...
    (2 + (4 + (AR * beta / eta) ^ 2 * ...
    (1 + (tand(sweep) / beta) ^ 2)) ^ 0.5) / Sref;

    CLa_supersonic = 4 * AR * Sexp * F /...
    (2 + (4 - (AR * beta / eta) ^ 2 * ...
    (1 + (tand(sweep) / beta) ^ 2)) ^ 0.5) / Sref;
   
    if M < 1
        CLa = [CLa CLa_subsonic];
        M_list = [M_list M];
    elseif M > 1 &&  M < 1 / cosd(le_sweep)
        
    else
        CLa = [CLa, CLa_supersonic];
    end
end

M = 0:0.01:2;
plot(M, CLa);
ylim([0 7])