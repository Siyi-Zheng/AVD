clear
clc

moment_x = 16e+06; % bending moment about x (Nm)
moment_y = 10; % bending moment about y (Nm)
d_fuslg = 6.3; % fuselage diameter (m)
circum_fuslg = pi*d_fuslg; % fuselage circumference (m)
no_stringer = 82; % number of stringers
dl_stringer = circum_fuslg/no_stringer; % stringer spacing
fprintf("Stringer Spacing : %.3g cm, %.3g inches \n",dl_stringer*100,dl_stringer*39.3701)
t_skin = 3e-3; % skin thickness (m)
A_stringer = 50e-6; % one stringer area (m^2)
dA_boom = ones(1,no_stringer) * 15*t_skin*t_skin; % skin collaborative area
A_boom =  A_stringer + dA_boom; % boom area (m^2)
fprintf("Boom Area : %.3g , Stringer Area: %.3g\n",A_boom(1),A_stringer);


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

% % plot boom location
% figure(1)
% plot(x,y,".",MarkerSize=8)
% 
% % Plot fuselage cross-section
% figure(2);
% hold on;
% grid on;
% plot3(x, y, z, 'LineWidth', 2,'Color','r'); % Circular fuselage cross-section
% plot3(x, y, z, '.','MarkerSize', 12,'Color','b') 
% % Plot stress vectors at each stringer location
% quiver3(x, y, z, u, v, sigx+sigy, 'k', 'LineWidth', 1, 'MaxHeadSize', 0.5);
% % Labels and title
% xlabel('x');
% ylabel('y');
% zlabel('Direct stress \sigma_z (MPa)');
% title('Direct stress distribution around fuselage');
% legend({'Fuselage cross-section', 'String Location','Direct Stress'}, 'Location', 'northeast');
% view(3); 

% using Z shape stringer
syms t_stringer; % stringer thickness
syms w_flange;  % stringer flange width
syms h_web; % stringer web height
Area = 2*w_flange*t_stringer + (h_web - 2*t_stringer)*t_stringer == A_stringer; % total area of stringer, need to be the same as A_stringer
I = (t_stringer*h_web^3)/12 + 2*(t_stringer*w_flange*(h_web/2)^2); % second moment of area of stringer


f = solve(Area, h_web); % Solve for h_web in terms of t_stringer and w_flange
f = matlabFunction(f, 'Vars', [t_stringer, w_flange]); % Convert symbolic function to numerical function
figure(3)
fsurf(f, [1e-3 5e-3 5e-2 20e-2]); % Plot h_web as a surface against t_stringer (1-4mm) and w_flange (5-20cm)
xlabel("Thickness (m)")
ylabel("Flange Width (m)")
zlabel("Web Height (m)")
title("Stringer Total Area")

t_fixed = 1e-3; % Fixed t_stringer value
I_func = subs(I, t_stringer, t_fixed);
I_func = matlabFunction(I_func, 'Vars', [w_flange, h_web]);
% % Plot the surface
figure(4);
fsurf(I_func, [1e-2 10e-2 5e-2 20e-2]);
xlabel('h_{web} (Height of Web)');
ylabel('w_{flange} (Width of Flange)');
zlabel('I (Moment of Inertia)');
title(sprintf('Moment of Inertia with t_{stringer} = %.1g mm', t_fixed*1000));
colorbar;
colormap(hot(20));
shading interp;

%%

% Define the objective function (negative for maximization)
objective = @(x) -((x(1)*x(2)^3)/12 + 2*(x(1)*x(3)*(x(2)/2)^2));

% Define the constraint function (equality constraint)
constraint = @(x) deal([], 2*x(3)*x(1) + (x(2) - 2*x(1))*x(1) - A_stringer);

% Set initial guesses [t_stringer, h_web, w_flange]
x0 = [2e-3, 0.1, 0.06];

% Define bounds (optional but helps convergence)
lb = [1e-3, 0.05, 0.02]; % Lower bounds
ub = [5e-3, 0.20, 0.15]; % Upper bounds

% Run optimization
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter');
[x_opt, I_max] = fmincon(objective, x0, [], [], [], [], lb, ub, constraint, options);

% Display results
disp('Optimal t_stringer, h_web, w_flange:');
disp(x_opt);
disp('Maximum I:');
disp(-I_max); % Since we minimized -I, negate it to get max value

%%



% check for failure
% material: Al 2024-T861
sig_yield = 431e6; % yield stress
fprintf("Maximum stress: %.2g MPa, Yield stress: %.2g MPa \n",sig_max/1e6,sig_yield/1e6);

% check for stringer buckling
E = 73.85e9; % Young's Modulus
L_eff = [0.8, 1.6, 3.2, 4.8, 6.4]; % effective length 
sig_euler = pi^2*E*(-I_max)./(L_eff.^2*A_stringer); % euler buckling stress
% sig_euler = double(subs(sig_euler, [t_stringer, w_flange, h_web, L_eff], [1.5e-3, 10e-3, 9e-3, 0.4]));
fprintf("Stringer buckling stress: %.2g MPa, %.2g MPa, %.2g MPa, %.2g MPa, %.2g MPa \n" ,sig_euler/1e6);

% check for skin buckling
sig_skin = 3.6*E*(t_skin/dl_stringer)^2;
fprintf("Skin buckling stress: %.2g MPa \n" ,sig_skin/1e6);



