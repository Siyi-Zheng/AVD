meshSizes = [0.82 1 2 3]; % mm
meshDensities = 1./meshSizes; % mm^-1
stresses = [240 236 213 196]; % MPa
eigs = [2.404 2.38 2.20 1.90]; % eigs

% extrapolated data created with the curve fitter tool
X = linspace(1/3, 4, 100);
stress2 = 243.44 ./ (1 + exp(-3.115 * (X + 0.1226)));
eigs2 = (-0.0081 * X.^2 + 2.5135 * X - 0.5131) ./ (X - 0.1629);

% plotting the data
ax = gca();
ha = plot(meshDensities, stresses, "kx-", LineWidth=2);
hold on
plot(X, stress2, "k--")
ylabel("Maximum stress (MPa)")
yyaxis right
hb = plot(meshDensities, eigs, "ro-", LineWidth=2)
plot(X, eigs2, "r--")
ylabel("Buckling eigenvalue")
xlabel("Mesh element density (1/mm)")
grid on
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = [0 0 0];
xlim([0.3 3])
legend([ha, hb], "Max. stress", "Buckling eigenvalues")
