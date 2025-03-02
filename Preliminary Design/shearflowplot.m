function q = shearflowplot(P , Q , T)



%data for fuselage ring
r = 3.17;


phi = 0:10:360;
phi = phi .* (pi/180);
rho = zeros(length(phi) , 1);
rho(:) = r;


%Calculating shear flow around the fuselage ring
q = (T + P*r)./(2 * pi * r^2) + (P .* cos(phi)) ./ (pi * r) + (Q .* sin(phi)) ./ (pi * r);


end