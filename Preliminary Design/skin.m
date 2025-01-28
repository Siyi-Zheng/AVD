% if there is an error run wing_load_distribution.m first
momentMax = max(M_x); % Nm
momentMin = min(M_x); % Nm
c = (0.7 - 0.12) .* chord; % m
b2 = 0.1255 .* chord; % m
E = 70; % GPa

% max. bending moment - no stringers
N = momentMax / (c .* b2); % N/m
t2 = (N ./ (3.62 * E ./ c .^ 2)) .^ (1/3); % mm
sigma = N ./ (t2 * 1000); % MPa
panelCSA = c .* t2 * 1000; % mm^2

% max bending moment - with stringers
n = 20;
b = c / n * 1000; % mm
ts = 1; % mm
h = 20; % mm
flangeWebRatio = 0.3; % ratio of flange length to web height
d = h * flangeWebRatio; % mm
As = ts * (h + 2 * d); % mm^2
t2 = (N ./ (3.62 * E ./ ...
    (b / 1000) .^ 2)) .^ (1/3); % mm
sigma_0 = N ./ (t2 * 1000); % MPa
t_e = t2 + (As ./ b); % mm
effectivePanelCSA = n .* b .* t_e; % mm^2

% catchpole diagram
catchpoleX = ts ./ t2; % x-value to read
catchpoleY = As ./ (b .* t2); % y-value to read
stressRatio = 1; % READ FROM THE DIAGRAM!!!