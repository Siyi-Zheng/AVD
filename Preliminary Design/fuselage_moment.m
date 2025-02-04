clear
clc

moment_x = 1e5; % bending moment about x
moment_y = 2e5; % bending moment about y
d_fuslg = 6.3; % fuselage diameter
dl_stringer = 0.15; % stringer spacing
no_stringer = 70; % number of stringers
t_skin = 0.003; % skin thickness
A_stringer = 0.00055; % one stringer area
A_boom = A_stringer + 15*t_skin*t_skin; % boom area

boom_angle = linspace(0,360,no_stringer);
x = d_fuslg/2 * cos(deg2rad(boom_angle));
y = d_fuslg/2 * sin(deg2rad(boom_angle));
z = zeros(1,no_stringer);
u = zeros(1,no_stringer);
v = zeros(1,no_stringer);

figure(1)
plot(x,y,".",MarkerSize=8)

Ix = sum(A_boom*(y.*y)); % second moment of area about x
Iy = sum(A_boom*(x.*x)); % second moment of area about y

x_max = max(x); % max x distance from centre
y_max = max(y); % max y distance from centre

sigx = moment_x*y/Ix;
sigy = moment_y*x/Iy;


sigx_max = max(sigx);
sigy_max = max(sigy);

% Plot fuselage cross-section
figure;
hold on;
grid on;
plot3(x, y, z, 'LineWidth', 2,'Color','r'); % Circular fuselage cross-section
plot3(x, y, z, '.','MarkerSize', 12,'Color','b') 

% Plot stress vectors at each stringer location
quiver3(x, y, z, u, v, sigx+sigy, 'k', 'LineWidth', 1, 'MaxHeadSize', 0.5);

% Labels and title
xlabel('x');
ylabel('y');
zlabel('Direct stress \sigma_z (MPa)');
title('Direct stress distribution around fuselage');
legend({'Fuselage cross-section', 'String Location','Direct Stress'}, 'Location', 'northeast');

view(3); % Set 3D view
hold off;

