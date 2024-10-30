clear all
clc

% read the data from polar_data.csv
% columns are alpha(deg), cl, cd, cdp, cm, top_xtr, bot_xtr, then repeats for the second airfoil
% remove the header row
% airfoils are in order NASA SC(2)-0714 (landing), NASA SC(2)-0614 (landing), NASA SC(2)-0714 (cruise), NASA SC(2)-0614 (cruise)
data = csvread('polar_data.csv', 1, 0);
data = data(1:end, :);
alpha = data(:, 1);
cl = data(:, 2);
cd = data(:, 3);
cdp = data(:, 4);
cm = data(:, 5);
top_xtr = data(:, 6);
bot_xtr = data(:, 7);
alpha2 = data(:, 8);
cl2 = data(:, 9);
cd2 = data(:, 10);
cdp2 = data(:, 11);
cm2 = data(:, 12);
top_xtr2 = data(:, 13);
bot_xtr2 = data(:, 14);
alpha3 = data(:, 15);
cl3 = data(:, 16);
cd3 = data(:, 17);
cdp3 = data(:, 18);
cm3 = data(:, 19);
top_xtr3 = data(:, 20);
bot_xtr3 = data(:, 21);
alpha4 = data(:, 22);
cl4 = data(:, 23);
cd4 = data(:, 24);
cdp4 = data(:, 25);
cm4 = data(:, 26);
top_xtr4 = data(:, 27);
bot_xtr4 = data(:, 28);

% prandtl-glauert correction at cruise (aerfoil 3 and 4)
M = 0.83;
beta = sqrt(1 - M^2);
cl3 = cl3 ./ beta;
cd3 = cd3 ./ beta;
cl4 = cl4 ./ beta;
cd4 = cd4 ./ beta;

% plot the cl vs alpha
figure(1)
plot(alpha, cl, 'k', 'LineWidth', 2)
hold on
plot(alpha2, cl2, 'r', 'LineWidth', 2)
plot(alpha3, cl3, 'r:', 'LineWidth', 2)
plot(alpha4, cl4, 'k:', 'LineWidth', 2)
xlabel('Angle of attack (deg)')
ylabel('Lift coefficient')
grid on
set(gca, 'FontSize', 12)
legend(["NASA SC(2)-0714 (landing)", "NASA SC(2)-0614 (landing)", "NASA SC(2)-0714 (cruise)", "NASA SC(2)-0614 (cruise)"], "Location", "northwest")

% plot the cd vs alpha
figure(2)
plot(alpha, cd, 'k', 'LineWidth', 2)
hold on
plot(alpha2, cd2, 'r', 'LineWidth', 2)
plot(alpha3, cd3, 'r:', 'LineWidth', 2)
plot(alpha4, cd4, 'k:', 'LineWidth', 2)
xlabel('Angle of attack (deg)')
ylabel('Drag coefficient')
grid on
set(gca, 'FontSize', 12)
legend(["NASA SC(2)-0714 (landing)", "NASA SC(2)-0614 (landing)", "NASA SC(2)-0714 (cruise)", "NASA SC(2)-0614 (cruise)"], "Location", "northwest")

% plot cl/cd vs alpha
figure(3)
plot(alpha, cl./cd, 'k', 'LineWidth', 2)
hold on
plot(alpha2, cl2./cd2, 'r', 'LineWidth', 2)
plot(alpha3, cl3./cd3, 'r:', 'LineWidth', 2)
plot(alpha4, cl4./cd4, 'k:', 'LineWidth', 2)
xlabel('Angle of attack (deg)')
ylabel('Lift-to-drag ratio')
grid on
set(gca, 'FontSize', 12)
legend(["NASA SC(2)-0714 (landing)", "NASA SC(2)-0614 (landing)", "NASA SC(2)-0714 (cruise)", "NASA SC(2)-0614 (cruise)"], "Location", "northwest")
