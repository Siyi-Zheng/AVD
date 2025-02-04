clear
clc



%%
%skin thickness shear flow
%need to break load cases down into 



%data for fuselage ring
r = 6.34;


phi = 0:10:360;
phi = phi .* (pi/180);
rho = zeros(length(phi) , 1);
rho(:) = r;

%tangiential load
P = 0; %N

%Radial Loading
Q = 2500000; %N

%Pure Torque
T = 0; %Nm

%Calculating shear flow around the fuselage ring
q = (T + P*r)/(2 * pi * r^2) + (P .* cos(phi)) ./ (pi * r) + (Q .* sin(phi)) ./ (pi * r);

%make positive, normalise and add the radius for visual appeal and ease of reading
q(:) = abs(q(:)/max(q)) + r;

%fuselage ring
fus_ring = zeros(length(phi) , 1);
fus_ring(:) = min(q);

%plotting shear flow distribution
figure

polarplot(phi , rho , Color='r' , LineWidth=2)
hold on
polarplot(phi , q , Color='c' , LineWidth=2)

ax = gca;
d = ax.ThetaDir;
ax.ThetaDir = 'counterclockwise';
ax.ThetaZeroLocation = 'bottom';

hold off


%%
%Pressurization Loads
%taking pressure difference at cruise
Ho = 43000;
Hi = 8000;
Ho = convlength(Ho,'ft','m');
Hi = convlength(Hi,'ft','m');
format long g;
[~,~,Po] = atmosisa(Ho);
[~,~,Pi] = atmosisa(Hi);

P = Pi - Po;

%fuselage is perfect cylinder with diameter 
D = 6.34;

stress_thickness_ratio_h = (P * D) / 2;

stress_thickness_ratio_l = (P * D) / 4;



