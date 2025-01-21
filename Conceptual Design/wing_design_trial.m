clear
clc
close all

m= (0.9-0.74)/(4-15.1);
airfoil_M_crit= @(t_c) ((0.9-0.74)/(4-15.1))*(t_c - 4) + 0.9;
fplot(airfoil_M_crit, [4 15.1])

airfoil_M_crit(14);
airfoil_M_crit(5);

