meshSizes = [0.75 1 2 3]; % mm
meshDensities = 1./meshSizes; % mm^-1
stresses = [229 222 201 194]; % MPa
eigs = [1.43 1.41 1.31 0.747]; % eigs

% plotting the data
plot(meshDensities, stresses, "kx-");
hold on
ylabel("Maximum stress (MPa)")
yyaxis right
plot(meshDensities, eigs, "rx-")
ylabel("Buckling eigenvalue")
xlabel("Mesh element density (1/mm)")
grid on
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = [0 0 0];

% extrapolated data created with the curve fitter tool
