clear
clc

table = readtable("catchpole_data.csv");
catchpoledata = table(1:8, 1:10);
farrardata= table(9:16, 1:10);

% % max. bending moment - no stringers

% t2 = (N ./ (3.62 * E ./ c .^ 2)) .^ (1/3) *1000; % mm
% sigma = N ./ (t2 * 1000); % MPa Critical buckling stress
% panelCSA = c .* t2 * 1000; % mm^2


%%%%%%%%%Design panel size and stringer size
% max bending moment - with stringers

%fixed variables
% if there is an error run wing_load_distribution.m first
momentMax = 4362716978.29012/100; % Nm
c = (0.7 - 0.12) .* 12.05; % m
b2 = 1.51; % m
E = 70; % GPa
N = momentMax / (c .* b2); % N/m compressive load per unit length
flangeWebRatio = 0.3; % ratio of flange length to web height

n = 20: 1: 50;
% ts = linspace(1,15, length(n)); % mm      %stringer thickness
% h = linspace(20, 300, length(n)); % mm      %stringer web height
ts= 1:0.1:15;
h= 100:5:300;
b = c ./ n.* 1000; % mm
d = h .* flangeWebRatio; % mm  %flange width of stringer

for i = 1:length(n)
    for j=1:length(ts)
        for k = 1:length(h)
            As(i,j,k)= ts(j) .* (h(k) + 2 * d(k)); % mm^2
            t2(i,j,k)= (N /3.62/ (E*10^9) * (b(i) / 1000) ^ 2) ^ (1/3) *1000; % mm
            sigma_0(i,j,k) = N / (t2(i,j,k) * 1000); % MPa
            t_e(i,j,k) = t2(i,j,k) + (As(i,j,k) / b(i)); % mm
            effectivePanelCSA = n(i)* b(i)* t_e(i,j,k); % mm^2
            % catchpole diagram
            catchpoleY(i,j,k) = ts(j)/ t2(i,j,k); % x-value to read
            catchpoleX(i,j,k) = As(i,j,k)/ (b(i)* t2(i,j,k)); % y-value to read
        end
    end
end

catchpoleY(catchpoleY > 1.4 | catchpoleY < 0.4) = 0;
catchpoleX(catchpoleX>1.2) = 0;
Y_grid = [ 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4]; % 1x10 (Y-coordinates)
X_grid = [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2];  % 8x1 (X-coordinates), NaN ignored
catchpoledata_matrix = table2array(catchpoledata);
farrardata_matrix = table2array(farrardata);
stressFactors = catchpoledata_matrix(2:end,2:end);
f_testing = farrardata_matrix(2:end,2:end);

for i = 1:length(n)
    for j=1:length(ts)
        for k = 1:length(h)
            if catchpoleY(i,j,k) ==0 || catchpoleX(i,j,k) ==0
                continue
            end
            stressRatio(i,j,k) = interp2(Y_grid,X_grid, stressFactors, catchpoleX(i,j,k),catchpoleY(i,j,k)); 
            F_factor(i,j,k) = interp2(Y_grid,X_grid, f_testing, catchpoleX(i,j,k),catchpoleY(i,j,k));
        end
    end
end

sigma_cr = stressRatio .*sigma_0;

%Farrar efficiency (stringer local buckling)

Rib_spacing = ; % m
for g = 1:length(Rib_spacing)
                    sigma_f(i,j,k)= F_factor(i,j,k)/sqrt(Rib_spacing(g)/N/(E*10^9))/10^6;

resultant = sigma_cr - sigma_f; 



