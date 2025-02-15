clear
clc

moment_x = 2.6393e+09; % bending moment about x (Nm)
moment_y = 1000; % bending moment about y (Nm)
d_fuslg = 6.3; % fuselage diameter (m)
circum_fuslg = pi*d_fuslg; % fuselage circumference (m)
no_stringer = 80; % number of stringers
dl_stringer = circum_fuslg/no_stringer; % stringer spacing
fprintf("Stringer Spacing : %.3g inches \n",dl_stringer*39.3701)
t_skin = 0.003; % skin thickness (m)
A_stringer = 60e-6; % one stringer area (m^2)
dA_boom = ones(1,no_stringer) * 15*t_skin*t_skin; % skin collaborative area
A_boom =  A_stringer + dA_boom; % boom area (m^2)
fprintf("Boom Area (0) : %.3g \n",A_boom(1));


boom_angle = linspace(0,360,no_stringer+1);
boom_angle = boom_angle(1:no_stringer);
x = d_fuslg/2 * cos(deg2rad(boom_angle)); % x location of booms
y = d_fuslg/2 * sin(deg2rad(boom_angle)); % y location of booms
z = zeros(1,no_stringer);
u = zeros(1,no_stringer);
v = zeros(1,no_stringer);

Ix = sum(A_boom.*(y.*y)); % second moment of area about x
Iy = sum(A_boom.*(x.*x)); % second moment of area about y

x_max = max(x); % max x distance from centre
y_max = max(y); % max y distance from centre

sigx = moment_x*y/Ix; % direct stress due to x bending moment
sigy = moment_y*x/Iy; % direct stress due to y bending moment
sig = sigx + sigy; % direct stress due to x and y bending moment


dA_boom(1) = t_skin/6 * (2 + sig(no_stringer)/sig(1)) + t_skin/6 * (2 + sig(2)/sig(1));
for i = 2:no_stringer-1
    dA_boom(i) = t_skin/6 * (2 + sig(i-1)/sig(i)) + t_skin/6 * (2 + sig(i+1)/sig(i));
end
dA_boom(no_stringer) = t_skin/6 * (2 + sig(no_stringer-1)/sig(no_stringer)) + t_skin/6 * (2 + sig(1)/sig(no_stringer));
A_boom =  A_stringer + dA_boom;
fprintf("Boom Area (1) : %.3g \n",A_boom(2));

for i = 1:5
    Ix = sum(A_boom.*(y.*y)); % second moment of area about x
    Iy = sum(A_boom.*(x.*x)); % second moment of area about y
    sigx = moment_x*y/Ix; % direct stress due to x bending moment
    sigy = moment_y*x/Iy; % direct stress due to y bending moment
    sig = sigx + sigy; % direct stress due to x and y bending moment
    dA_boom(1) = t_skin/6 * (2 + sig(no_stringer)/sig(1)) + t_skin/6 * (2 + sig(2)/sig(1));
    for j = 2:no_stringer-1
        dA_boom(j) = t_skin/6 * (2 + sig(j-1)/sig(j)) + t_skin/6 * (2 + sig(j+1)/sig(j));
    end
    dA_boom(no_stringer) = t_skin/6 * (2 + sig(no_stringer-1)/sig(no_stringer)) + t_skin/6 * (2 + sig(1)/sig(no_stringer));
    A_boom =  A_stringer + dA_boom;
    fprintf("Boom Area (%d) : %.3g \n",i+1,A_boom(2));
end


sig_max = max(sig);

% plot boom location
figure(1)
plot(x,y,".",MarkerSize=8)

% Plot fuselage cross-section
figure(2);
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


% using Z shape stringer
syms t_stringer; % stringer thickness
syms w_flange;  % stringer flange width
syms h_web; % stringer web height
Area = 2*w_flange*t_stringer + (h_web - 2*t_stringer)*t_stringer == A_stringer; % total area of stringer
% need to be the same as A_stringer
I = (t_stringer*h_web^3)/12 + 2*(t_stringer*w_flange*(h_web/2)^2); % second moment of area of stringer
f = solve(Area,h_web);
f = matlabFunction(f, 'Vars', [t_stringer, w_flange]);
fsurf(f,[1e-3 5e-3 5e-2 20e-2])


% check for failure
% material: Al 2024-T861
sig_yield = 431e6; % yield stress
fprintf("Maximum stress: %.2g MPa \nYield stress: %.2g MPa",sig_max/1e6,sig_yield/1e6);

% check for stringer buckling
E = 73.85e9; % Young's Modulus
syms L_eff % effective length
sig_euler = pi^2*E*I/(L_eff^2*A_stringer); % euler buckling stress

% check for skin buckling
sig_skin = 3.6*E*(t_skin/dl_stringer)^2;


