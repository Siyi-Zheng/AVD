function [sig_skin, sig_euler, sig_max, Af_min] = fuslg_failure(no_stringer,A_stringer,Lfs)

max_bending = 15975543.3205404;
t_skin = 2e-3; % skin thickness (m)
moment_x = max_bending; % bending moment about x (Nm)
moment_y = 10; % bending moment about y (Nm)
d_fuslg = 6.3; % fuselage diameter (m)
circum_fuslg = pi*d_fuslg; % fuselage circumference (m)
dl_stringer = circum_fuslg/no_stringer; % stringer spacing
dA_boom = ones(1,no_stringer) * 15*t_skin*t_skin; % skin collaborative area
A_boom =  A_stringer + dA_boom; % boom area (m^2)


boom_angle = linspace(0,360,no_stringer+1);
boom_angle = boom_angle(1:no_stringer);
x = d_fuslg/2 * cos(deg2rad(boom_angle)); % x location of booms
y = d_fuslg/2 * sin(deg2rad(boom_angle)); % y location of booms

Ix = sum(A_boom.*(y.*y)); % second moment of area about x
Iy = sum(A_boom.*(x.*x)); % second moment of area about y


sigx = moment_x*y/Ix; % direct stress due to x bending moment
sigy = moment_y*x/Iy; % direct stress due to y bending moment
sig = sigx + sigy; % direct stress due to x and y bending moment


dA_boom(1) = t_skin/6 * (2 + sig(no_stringer)/sig(1)) + t_skin/6 * (2 + sig(2)/sig(1));
for i = 2:no_stringer-1
    dA_boom(i) = t_skin/6 * (2 + sig(i-1)/sig(i)) + t_skin/6 * (2 + sig(i+1)/sig(i));
end
dA_boom(no_stringer) = t_skin/6 * (2 + sig(no_stringer-1)/sig(no_stringer)) + t_skin/6 * (2 + sig(1)/sig(no_stringer));
A_boom =  A_stringer + dA_boom;

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
end

sig_max = max(sig);


% using Z shape stringer
% Define the objective function (negative for maximization)
objective = @(x) -((x(1)*x(2)^3)/12 + 2*(x(1)*x(3)*(x(2)/2)^2));

% Define the constraint function (equality constraint)
constraint = @(x) deal([], 2*x(3)*x(1) + (x(2) - 2*x(1))*x(1) - A_stringer);

% Set initial guesses [t_stringer, h_web, w_flange]
x0 = [2e-3, 0.1, 0.06];

% Define bounds [t_stringer, h_web, w_flange]
lb = [1e-3, 0.02, 0.02]; % Lower bounds (A = 58e-6)
ub = [5e-3, 0.10, 0.10]; % Upper bounds (A = 1450e-6)

% Run optimization
options = optimoptions('fmincon', 'Algorithm', 'sqp');
[x_opt, I_max] = fmincon(objective, x0, [], [], [], [], lb, ub, constraint, options);



D = 6.34; %frame diameter
Cf = 1/16000; % empirical based off work performed by Lockheed Martin and Shanley, J. Aero Science
M_ult = 1.6e7; % max bending stress
E = 73.85e9; % Young's Modulus 


%set a fixed IF_xx
If = (Cf * M_ult * D^2) / (E * Lfs);


% Define the objective function (minimizing Af)
objective = @(x) (2 * x(2) + x(1)) * x(3);

% Define the constraint function (equality constraint)
constraint = @(x) deal([], x(3) * ((x(1)^3 / 12) + (0.5 * x(2) * x(1)^2)) - If);

% Set initial guesses [h, b, t]
x0 = [0.15, 0.03, 2e-3];

% Define bounds [h, b, t]
lb = [0.10, 0.02, 1e-3]; % Lower bounds (I = 1.833e-7, 2.97m)
ub = [0.20, 0.1, 5e-3]; % Upper bounds  (I = 1.333e-5, 0.04m)

% Run optimization
options = optimoptions('fmincon', 'Algorithm', 'sqp');
[x_opt, Af_min] = fmincon(objective, x0, [], [], [], [], lb, ub, constraint, options);


% check for stringer buckling
sig_euler = pi^2*E*(-I_max)/(Lfs^2*A_stringer); % euler buckling stress


% check for skin buckling
w = 1.90*t_skin*sqrt(E/sig_max);
if w < dl_stringer*0.8
    sig_skin = 3.6*E*(t_skin/(dl_stringer - w))^2;
else
    sig_skin = 3.6*E*(t_skin/(dl_stringer))^2;
end



end