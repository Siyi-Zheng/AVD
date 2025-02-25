angles = [0 1 2 3]; % deg
stresses = [240 249 266 326]; % MPa, 1.126

% plotting the data
plot(angles, stresses, "kx", LineWidth=2);
xlabel("Load angle (deg)")
hold on
ylabel("Maximum stress (MPa)")
grid on
xlim([0 3.1])

% curve fitter
X = linspace(0, 3.1, 100);
Y = 239.66 * exp(0.0299*X) + 0.3845 * exp(1.7044*X);
plot(X, Y, "k--")
yline(276)