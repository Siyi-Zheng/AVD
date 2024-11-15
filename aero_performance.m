clear
clc

AR = 8.77;          %aspect ratio
Cla = 7.16;         %lift coeff. of aerofoil
d = 6.34;           %diameter of fuselage
b = 65;             %wing span
Sexp = 440;         %exposed area
Sref = 482;         %reference area
sweep = 24;         % sweep at max thickness pt.
le_sweep = 30;      %leading edge

CLa1 = [];
CLa2 = [];
CLa3 = [];
M_list1 = [];
M_list2 = [];
for M = 0:0.01:2

    beta = (1 - M ^ 2) ^ 0.5;
    eta = beta * Cla / (2 * pi);
    F = 1.07 * (1 + d / b) ^ 2;

    CLa_subsonic3 = 2 * pi * AR * Sexp * F /(2 + (4 + (AR * beta / eta) ^ 2 *(1 + (tand(sweep) / beta) ^ 2)) ^ 0.5) / Sref;

    %CLa_subsonic = (2*pi*AR*F*(Sexp/Sref)) / (2+(sqrt(4 + ((AR*beta/eta)^2)*(1+((tand(sweep))^2)/(beta*beta)))));
    CLa_subsonic = 2*pi/sqrt(1-M^2);

    CLa_supersonic = 4 * AR * Sexp * F /(2 + (4 + (AR * beta / eta) ^ 2 * (1 + (tand(sweep) / beta) ^ 2)) ^ 0.5) / Sref;
   %CLa_supersonic = 4/sqrt(-1+M^2);

    if M < 1
        CLa1 = [CLa1 CLa_subsonic];
        CLa3 = [CLa3 CLa_subsonic3];
        M_list1 = [M_list1 M];
    elseif M > 1
        CLa2 = [CLa2, CLa_supersonic];
        M_list2 = [M_list2 M];
    end
end

figure
plot(M_list1, CLa1,"r");
hold on
plot(M_list2, CLa2,"b");
plot(M_list1, CLa3,"g");
plot([0.85 0.85],[0,10],"k");
plot([1.2 1.2],[0,10],"k");
plot([1 1],[0,10],"k");
xlabel("M");
ylabel("CLa")
ylim([0 10]);
legend("Subsonic","Supersonic");
grid on









%IMPORTANT VALUES!!!!!
%CLa = 5.38124 FOR THE WING FOR M = 0.83