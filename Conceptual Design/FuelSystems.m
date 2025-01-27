clc
clear
close all

%Wing available volume calcualtor

%Need to get the available volume as a fuction of chord length

%airfoil data
 airfoil =  [0.000000,  0.000000;
  0.002000 , 0.010800;
  0.005000 , 0.016600;
  0.010000 , 0.022500;
  0.020000 , 0.029800;
  0.030000 , 0.034900;
  0.040000 , 0.038700;
  0.050000 , 0.041800;
  0.060000 , 0.044500;
  0.070000 , 0.046800;
  0.080000 , 0.048900;
  0.090000 , 0.050800;
  0.100000 , 0.052500;
  0.110000 , 0.054100;
  0.120000 , 0.055500;
  0.130000 , 0.056800;
  0.140000 , 0.058000;
  0.150000 , 0.059100;
  0.160000 , 0.060200;
  0.170000 , 0.061200;
  0.180000 , 0.062100;
  0.190000 , 0.062900;
  0.200000 , 0.063700;
  0.210000 , 0.064400;
  0.220000 , 0.065100;
  0.230000 , 0.065700;
  0.240000 , 0.066300;
  0.250000 , 0.066800;
  0.260000 , 0.067300;
  0.270000 , 0.067800;
  0.280000 , 0.068200;
  0.290000 , 0.068600;
  0.300000 , 0.068900;
  0.310000 , 0.069200;
  0.320000 , 0.069400;
  0.330000 , 0.069600;
  0.340000 , 0.069800;
  0.350000 , 0.069900;
  0.360000 , 0.070000;
  0.370000 , 0.070100;
  0.380000 , 0.070100;
  0.390000 , 0.070100;
  0.400000 , 0.070100;
  0.410000 , 0.070000;
  0.420000 , 0.069900;
  0.430000 , 0.069800;
  0.440000 , 0.069600;
  0.450000 , 0.069400;
  0.460000 , 0.069200;
  0.470000 , 0.069000;
  0.480000 , 0.068700;
  0.490000 , 0.068400;
  0.500000 , 0.068100;
  0.510000 , 0.067700;
  0.520000 , 0.067300;
  0.530000 , 0.066900;
  0.540000 , 0.066400;
  0.550000 , 0.065900;
  0.560000 , 0.065300;
  0.570000 , 0.064700;
  0.580000 , 0.064000;
  0.590000 , 0.063300;
  0.600000 , 0.062600;
  0.610000 , 0.061800;
  0.620000 , 0.061000;
  0.630000 , 0.060100;
  0.640000 , 0.059100;
  0.650000 , 0.058100;
  0.660000 , 0.057000;
  0.670000 , 0.055900;
  0.680000 , 0.054700;
  0.690000 , 0.053500;
  0.700000 , 0.052200;
  0.710000 , 0.050900;
  0.720000 , 0.049500;
  0.730000 , 0.048100;
  0.740000 , 0.046600;
  0.750000 , 0.045100;
  0.760000 , 0.043600;
  0.770000 , 0.042000;
  0.780000 , 0.040400;
  0.790000 , 0.038700;
  0.800000 , 0.037000;
  0.810000 , 0.035200;
  0.820000 , 0.033400;
  0.830000 , 0.031600;
  0.840000 , 0.029700;
  0.850000 , 0.027800;
  0.860000 , 0.025800;
  0.870000 , 0.023800;
  0.880000 , 0.021800;
  0.890000 , 0.019700;
  0.900000 , 0.017600;
  0.910000 , 0.015400;
  0.920000 , 0.013200;
  0.930000 , 0.010900;
  0.940000 , 0.008600;
  0.950000 , 0.006200;
  0.960000 , 0.003800;
  0.970000 , 0.001300;
  0.980000 ,-0.001300;
  0.990000 ,-0.003900;
  1.000000 ,-0.006600;

  0.000000 , 0.000000;
  0.002000 ,-0.010800;
  0.005000 ,-0.016600;
  0.010000 ,-0.022500;
  0.020000 ,-0.029800;
  0.030000 ,-0.034900;
  0.040000 ,-0.038800;
  0.050000 ,-0.041900;
  0.060000 ,-0.044600;
  0.070000 ,-0.046900;
  0.080000 ,-0.049000;
  0.090000 ,-0.050900;
  0.100000 ,-0.052600;
  0.110000 ,-0.054200;
  0.120000 ,-0.055700;
  0.130000 ,-0.057000;
  0.140000 ,-0.058200;
  0.150000 ,-0.059400;
  0.160000 ,-0.060500;
  0.170000 ,-0.061500;
  0.180000 ,-0.062400;
  0.190000 ,-0.063300;
  0.200000 ,-0.064100;
  0.210000 ,-0.064800;
  0.220000 ,-0.065500;
  0.230000 ,-0.066100;
  0.240000 ,-0.066700;
  0.250000 ,-0.067200;
  0.260000 ,-0.067700;
  0.270000 ,-0.068100;
  0.280000 ,-0.068500;
  0.290000, -0.068800;
  0.300000, -0.069100;
  0.310000, -0.069300;
  0.320000, -0.069500;
  0.330000, -0.069700;
  0.340000, -0.069800;
  0.350000, -0.069900;
  0.360000, -0.069900;
  0.370000, -0.069800;
  0.380000, -0.069700;
  0.390000, -0.069600;
  0.400000, -0.069400;
  0.410000, -0.069200;
  0.420000, -0.068900;
  0.430000, -0.068600;
  0.440000, -0.068200;
  0.450000, -0.067700;
  0.460000, -0.067200;
  0.470000, -0.066600;
  0.480000, -0.065900;
  0.490000, -0.065100;
  0.500000, -0.064200;
  0.510000, -0.063200;
  0.520000, -0.062200;
  0.530000, -0.061100;
  0.540000, -0.059900;
  0.550000, -0.058600;
  0.560000, -0.057200;
  0.570000, -0.055700;
  0.580000, -0.054100;
  0.590000, -0.052500;
  0.600000, -0.050800;
  0.610000, -0.049100;
  0.620000, -0.047300;
  0.630000, -0.045500;
  0.640000, -0.043600;
  0.650000, -0.041700;
  0.660000, -0.039700;
  0.670000, -0.037700;
  0.680000, -0.035600;
  0.690000, -0.033600;
  0.700000, -0.031500;
  0.710000, -0.029400;
  0.720000, -0.027400;
  0.730000, -0.025300;
  0.740000, -0.023300;
  0.750000, -0.021300;
  0.760000, -0.019300;
  0.770000, -0.017400;
  0.780000, -0.015500;
  0.790000, -0.013700;
  0.800000, -0.011900;
  0.810000, -0.010200;
  0.820000, -0.008600;
  0.830000, -0.007200;
  0.840000, -0.005900;
  0.850000, -0.004700;
  0.860000, -0.003700;
  0.870000, -0.002900;
  0.880000, -0.002300;
  0.890000, -0.001900;
  0.900000, -0.001700;
  0.910000, -0.001700;
  0.920000, -0.001900;
  0.930000, -0.002400;
  0.940000, -0.003100;
  0.950000, -0.004100;
  0.960000, -0.005400;
  0.970000, -0.006900;
  0.980000, -0.008700;
  0.990000, -0.010800;
  1.000000, -0.006600];
% Separate upper and lower surfaces
upper_surfacea = airfoil(1:103, :);
lower_surfacea = airfoil(104:206, :);

% Compute the thickness (distance between upper and lower surfaces at each x position)
thickness = upper_surfacea(:, 2) - lower_surfacea(:, 2);

% Numerical integration to find area under the airfoil
area_under_airfoil = cumtrapz(upper_surfacea(:, 1), thickness);

% Total area of the airfoil with chord length of 1
final_area = area_under_airfoil(end);

%disp(["Final CSA for 1m of chord length is" num2str(final_area)])

% Define a range of chord lengths to calculate the CSA
chord_lengths = [0,1.0,1.5,2.0]; % Example chord lengths in meters

% Calculate CSA for each chord length and display
for i = 1:length(chord_lengths)
    chord = chord_lengths(i);
    CSA = final_area * chord^2;
    %disp(['For a chord length of ', num2str(chord), ' m, the cross-sectional area is: ', num2str(CSA), ' square meters']);
end

%calculating chord distribution

%leading edge sweep angle
b = 32.45 * tand(26.6);
a = b - 3.47/4;

aa = (a + 3.47)/32.45;

%sweep_angle = atand(aa);
sweep_angle = atand(aa);

%chord distribution for first section
outer_chord = 13.89 - (9.75 * tand(sweep_angle));

gradient1 = (outer_chord - 13.89)/9.75;

%for outer section
gradient2 = (3.47 - outer_chord)/(32.45 - 9.75);

span1 = 2.55:0.01:9.75;
span2 = 9.75:0.01:29.45;

c2 = outer_chord - gradient2 * 9.75;

chord1 = gradient1 * span1 + 13.89;
chord2 = gradient2 * span2 + c2;

%plotting chord and available csa distributions

figure
hold on
plot(span1, chord1)
plot(span2, chord2)
title("Chord Ratio")
hold off


%plotting area distribution
area1 = (final_area * chord1.^2) ./ cosd(sweep_angle);
area2 = (final_area * chord2.^2) ./ cosd(sweep_angle);

figure
hold on
plot(span1, area1)
plot(span2, area2)
title("Area Ratio")
hold off

%volume calculations
volume_available1 = (trapz(span1 , area1));
volume_available2 = trapz(span2 , area2);



total_availablevolume = volume_available1 + volume_available2;

%fuel volume required at 15 degrees density of 804kg/m^3

Wf = 1.760909890249436e+05;

rhof = 804;

volume_hold = Wf / rhof;

volume_holdl = volume_hold * 1000;

%increasing to factor foam area, tank area and extra volume for desnity

volume_req = volume_hold * 1.10;

%volume_wings = (volume_hold * 0.60) * 2;

% Fuel tank data (transposed for easier handling in rows)
%{
fueltank = [
    0.150000 , 0.059100;
    0.160000 , 0.060200;
    0.170000 , 0.061200;
    0.180000 , 0.062100;
    0.190000 , 0.062900;
    0.200000 , 0.063700;
    0.210000 , 0.064400;
    0.220000 , 0.065100;
    0.230000 , 0.065700;
    0.240000 , 0.066300;
    0.250000 , 0.066800;
    0.260000 , 0.067300;
    0.270000 , 0.067800;
    0.280000 , 0.068200;
    0.290000 , 0.068600;
    0.300000 , 0.068900;
    0.310000 , 0.069200;
    0.320000 , 0.069400;
    0.330000 , 0.069600;
    0.340000 , 0.069800;
    0.350000 , 0.069900;
    0.360000 , 0.070000;
    0.370000 , 0.070100;
    0.380000 , 0.070100;
    0.390000 , 0.070100;
    0.400000 , 0.070100;
    0.410000 , 0.070000;
    0.420000 , 0.069900;
    0.430000 , 0.069800;
    0.440000 , 0.069600;
    0.450000 , 0.069400;
    0.460000 , 0.069200;
    0.470000 , 0.069000;
    0.480000 , 0.068700;
    0.490000 , 0.068400;
    0.500000 , 0.068100;
    0.510000 , 0.067700;
    0.520000 , 0.067300;
    0.530000 , 0.066900;
    0.540000 , 0.066400;
    0.550000 , 0.065900;
    0.560000 , 0.065300;
    0.570000 , 0.064700;
    0.580000 , 0.064000;
    0.590000 , 0.063300;
    0.600000 , 0.062600;
    0.610000 , 0.061800;
    0.620000 , 0.061000;
    0.630000 , 0.060100;
    0.640000 , 0.059100;
    0.650000 , 0.058100;
    0.660000 , 0.057000;
    0.670000 , 0.055900;
    0.680000 , 0.054700;
    0.690000 , 0.053500;
    0.700000 , 0.052200;
    0.700000, -0.031500;
    0.690000, -0.033600;
    0.680000, -0.035600;
    0.670000, -0.037700;
    0.660000, -0.039700;
    0.650000, -0.041700;
    0.640000, -0.043600;
    0.630000, -0.045500;
    0.620000, -0.047300;
    0.610000, -0.049100;

    0.600000, -0.050800;
    0.590000, -0.052500;
    0.580000, -0.054100;
    0.570000, -0.055700;
    0.560000, -0.057200;
    0.550000, -0.058600;
    0.540000, -0.059900;
    0.530000, -0.061100;
    0.520000, -0.062200;
    0.510000, -0.063200;
    0.500000, -0.064200;
    0.490000, -0.065100;
    0.480000, -0.065900;
    0.470000, -0.066600;
    0.460000, -0.067200;
    0.450000, -0.067700;
    0.440000, -0.068200;
    0.430000, -0.068600;
    0.420000, -0.068900;
    0.410000, -0.069200;
    0.400000, -0.069400;
    0.390000, -0.069600;
    0.380000, -0.069700;
    0.370000, -0.069800;
    0.360000, -0.069900;
    0.350000, -0.069900;
    0.340000, -0.069800;
    0.330000, -0.069700;
    0.320000, -0.069500;
    0.310000, -0.069300;
    0.300000, -0.069100;
    0.290000, -0.068800;
    0.280000, -0.068500;
    0.270000, -0.068100;
    0.260000, -0.067700;
    0.250000, -0.067200;
    0.240000, -0.066700;
    0.230000, -0.066100;
    0.220000, -0.065500;
    0.210000, -0.064800;
    0.200000, -0.064100;
    0.190000, -0.063300;
    0.180000, -0.062400;
    0.170000, -0.061500;
    0.160000, -0.060500;
    0.150000, -0.059400
    0.150000 , 0.059100]';
%}

fueltank = [upper_surfacea(17:71 , :) ; flipud(lower_surfacea(17:71 , :)) ; upper_surfacea(17 , :)];
%central tank and front tank seperation
%structural member coordinates

% Separate upper and lower surfaces
upper_surface = fueltank(1:51 , :); % Columns 1-46 for upper surface
lower_surface = fueltank(52:102 , :); % Columns 47-92 for lower surface

% Calculate thickness between upper and lower surfaces
thickness = upper_surface(: , 2) - flipud(lower_surface(: , 2));

% Calculate the area by integrating thickness along the chord positions
Tankarea = trapz(upper_surface(: , 1), thickness);

% Display the area
%disp(['The total cross-sectional area of the fuel tank is: ', num2str(Tankarea), ' square units']);


%creating spar 
spar1 = [upper_surfacea(15:17,:) ; flipud(lower_surfacea(15:17,:)) ; upper_surfacea(15,:)];
spar2 = [upper_surfacea(71:73,:) ; flipud(lower_surfacea(71:73,:)) ; upper_surfacea(71,:)];



%plotting area distribution
area3 = (Tankarea * chord1.^2) / cosd(sweep_angle);
area4 = (Tankarea * chord2.^2) / cosd(sweep_angle);

figure
hold on
plot(span1, area3)
plot(span2, area4)
title("Area Ratio")
hold off

%volume calculations
volume1 = (trapz(span1 , area3) * 0.85);
volume2 = (trapz(span2 , area4) * 0.85);

total_tank_volume = (volume1 * 2 + volume2 * 2) - 13;

Volumer_needed = volume_req - total_tank_volume;

% Define the collector area by combining upper and lower surfaces
collector = [upper_surfacea(1:11, :); flipud(lower_surfacea(1:11, :)); upper_surfacea(1, :)];

% Separate collector into upper and lower parts
collector_upper = collector(1:11, :); % First 15 points for upper surface
collector_lower = collector(12:22, :); % Next 15 points for lower surface

% Calculate thickness between upper and lower surfaces (only y-coordinates)
thickness = collector_upper(:, 2) - flipud(collector_lower(:, 2));

% Calculate the area of the collector using numerical integration
collector_area = trapz(collector_upper(:, 1), thickness);

areac_ratio1 = collector_area * chord1.^2;
areac_ratio2 = collector_area * chord2.^2;

figure
hold on 
plot(span1 , areac_ratio1)
plot(span2 , areac_ratio2)
hold off

% Calculate collector area integration over specified ranges
collector1 = trapz(span1(240:676), areac_ratio1(240:676)); % Corrected typo
collector2 = trapz(span2(1195:1971), areac_ratio2(1195:1971));      % Corrected typo


%collector lengths
length1 = (span1(676) - span1(240)) * tand(30.15);
length2 = (span2(1971) - span2(1195)) * tand(30.15);


difference = Volumer_needed;


central_tank_length = difference / 8.165;

%%

x_slat = [0.0345,0.0325,0.0321,0.0343,0.039,0.0439,0.0675,0.075,0.1069];
y_slat = [-0.0368,-0.0283,-0.0185,-0.0008,0.0113,0.0225,0.0385,0.0425,0.0536];

%{
%slat calculation
% Given points
x = [0.1069, 0.0345, 0.0321]; % x-coordinates
y = [0.036, -0.0368, -0.0361]; % y-coordinates

% Construct the matrix system Ax = b
A = [sqrt(x(1)), x(1), 1;
     sqrt(x(2)), x(2), 1;
     sqrt(x(3)), x(3), 1];

b = [y(1); y(2); y(3)];

% Solve for coefficients [a; b; c]
coefficients = A\b;

% Extract coefficients
a = coefficients(1);
b = coefficients(2);
c = coefficients(3);

disp(['a = ', num2str(a), ', b = ', num2str(b), ', c = ', num2str(c)]);

% Generate x values for plotting
x_plot = linspace(min(x), max(x), 100);

% Compute the corresponding y values
y_plot = a * sqrt(x_plot) + b * x_plot + c;

% Define the endpoints of the line
x1 = 0.1069; y1 = 0.0536; % First endpoint
x2 = 0.0345; y2 = -0.0361; % Second endpoint

% Define the points to reflect
points = [x_plot' , y_plot']; % Example points (x, y) as rows

% Line equation coefficients
a = y2 - y1;
b = x1 - x2;
c = x2 * y1 - x1 * y2;

% Initialize matrix for reflected points
reflected_points = zeros(size(points));

% Reflect each point
for i = 1:size(points, 1)
    x = points(i, 1);
    y = points(i, 2);
    x_reflected = x - 2 * a * (a * x + b * y + c) / (a^2 + b^2);
    y_reflected = y - 2 * b * (a * x + b * y + c) / (a^2 + b^2);
    reflected_points(i, :) = [x_reflected, y_reflected];
end

x_slat = reflected_points(:,1);
y_slat = reflected_points(:,2);

%}

%%
%coordinates for flaps
x_flap = [0.2512,0.2501,0.2474,0.2424,0.2367,0.2297,0.2217,0.2146,0.209,0.2003,0.1892,0.1815,0.1759];
y_flap = [-0.021,-0.0108,-0.0015,0.0072,0.01320,0.0181,0.0219,0.0244,0.026,0.0279,0.0298,0.0307,0.0327];

x_flap = 1 - x_flap;
%%
%plotting

% Plot the airfoil
figure
figure('Position', [100, 100, 1000, 370]); % [x, y, width, height]
hold on

% Plot the upper and lower surfaces of the airfoil
plot(upper_surfacea(:, 1), upper_surfacea(:, 2), 'b', 'DisplayName', 'Upper Surface');
plot(lower_surfacea(:, 1), lower_surfacea(:, 2), 'r', 'DisplayName', 'Lower Surface');

% Plot the fuel tank outline
plot(fueltank(:,1), fueltank(:,2), 'k', 'HandleVisibility', 'off'); % Not in legend

% Fill the area of the fuel tank with transparent black
x_fill = fueltank(:,1);        % x-coordinates of the fuel tank section
y_fill = fueltank(:,2);        % y-coordinates of the fuel tank section
fill(x_fill, y_fill, 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Fuel Tank'); % Transparent black shading

% Plotting spars with gray color
plot(spar1(:,1), spar1(:,2), 'Color', [0.5, 0.5, 0.5], 'HandleVisibility', 'off'); % Not in legend
plot(spar2(:,1), spar2(:,2), 'Color', [0.5, 0.5, 0.5], 'HandleVisibility', 'off'); % Not in legend

% Fill area between spar1 and spar2 with transparent gray
xs_fill1 = spar1(:,1); % x-coordinates of spar1
ys_fill1 = spar1(:,2); % y-coordinates of spar1
xs_fill2 = spar2(:,1); % x-coordinates of spar2
ys_fill2 = spar2(:,2); % y-coordinates of spar2

% Combine coordinates to create a closed polygon for filling
x_fill_spar = [xs_fill1; flipud(xs_fill2)];
y_fill_spar = [ys_fill1; flipud(ys_fill2)];

% Fill the area between the spars
fill(x_fill_spar, y_fill_spar, [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Spar Fill');

%plot slats
plot(x_slat , y_slat, 'k-', 'LineWidth', 1 , 'DisplayName','Slat'); % Curve

%plotting flaps
plot(x_flap, y_flap , 'k-' , 'LineWidth', 1 , 'DisplayName','Flap') %curve for flap

% Additional plot settings
xlabel('Chord Position (x/c)');
ylabel('Thickness (y/c)');
%title('Wing Cross Section with Fuel Tank Integration');

% Maintain equal scaling and set axis limits
axis equal;
xlim([-0.1, 1.1]);
%ylim([-0.5, 0.5]);

legend;
axis equal
hold off

% Additional plot settings
xlabel('Chord Position (x/c)');
ylabel('Thickness (y/c)');
%title('Wing Cross Section with Fuel Tank Integration');
legend;
axis equal
hold off


%TAILPLANE TANK SIZING


% Airfoil data (with y-coordinates inverted)
lower_surfacet = [
    0.000000  -0.000000;
    0.000150  -0.009560;
    0.006280  -0.020300;
    0.018650  -0.031760;
    0.037300  -0.043240;
    0.062030  -0.053820;
    0.092300  -0.062650;
    0.127320  -0.069150;
    0.166040  -0.073200;
    0.207380  -0.075240;
    0.251310  -0.075970;
    0.297960  -0.075540;
    0.346810  -0.074020;
    0.397330  -0.071500;
    0.448970  -0.068110;
    0.501170  -0.063970;
    0.553350  -0.059240;
    0.604960  -0.054050;
    0.655410  -0.048540;
    0.704170  -0.042850;
    0.750700  -0.037120;
    0.794490  -0.031450;
    0.835060  -0.025970;
    0.871970  -0.020790;
    0.904820  -0.016020;
    0.933240  -0.011760;
    0.956930  -0.008120;
    0.975630  -0.005180;
    0.989140  -0.003020;
    0.997300  -0.001700;
    1.000030  -0.001260
];

upper_surfacet = [
    0.000000  -0.000000;
    0.005330  0.007920;
    0.015570  0.014010;
    0.030290  0.018700;
    0.049150  0.022480;
    0.071950  0.025860;
    0.098680  0.029220;
    0.129540  0.032820;
    0.164830  0.036600;
    0.204830  0.040160;
    0.248690  0.042830;
    0.295310  0.044460;
    0.344180  0.045100;
    0.394760  0.044820;
    0.446500  0.043710;
    0.498830  0.041880;
    0.551170  0.039450;
    0.602960  0.036550;
    0.653600  0.033270;
    0.702570  0.029750;
    0.749300  0.026070;
    0.793300  0.022350;
    0.834070  0.018660;
    0.871180  0.015120;
    0.904200  0.011800;
    0.932790  0.008800;
    0.956610  0.006210;
    0.975430  0.004100;
    0.989010  0.002540;
    0.997220  0.001580;
    0.999970  0.001260
];

% Plotting the inverted airfoil
figure;
hold on;

% Plot upper and lower surfaces with inverted y-coordinates
plot(upper_surfacet(:, 1), upper_surfacet(:, 2), 'b', 'LineWidth', 1, 'DisplayName', 'Upper Surface');
plot(lower_surfacet(:, 1), lower_surfacet(:, 2), 'r', 'LineWidth', 1, 'DisplayName', 'Lower Surface');

% Additional plot settings
xlabel('Chord Position (x/c)');
ylabel('Thickness (y/c)');
title('Inverted Airfoil Profile');
legend;
grid on;
axis equal;
hold off;

% parameters
area = 40.28 * 2; % m^2
AR = 5.8;
wing_span = (area * AR / 2) ^ 0.5; %m
taper_ratio = 0.4;
quarter_chord_sweep = 35; % degrees

% simple design
total_chord = 1 + 1 / taper_ratio;
tip_chord = 2 * area / (wing_span * total_chord);
root_chord = tip_chord / taper_ratio;
qcy = -wing_span / 2 * tand(quarter_chord_sweep);
le_tip = qcy + tip_chord / 4;
te_tip = qcy - 3 * tip_chord / 4;
le_root = root_chord / 4;
te_root = -3 * root_chord / 4;

chord_grad = (tip_chord - root_chord) / wing_span;

wing_span1 = 0:0.01:wing_span;

chordt = chord_grad * wing_span1 + root_chord;


% Calculate thickness (distance between upper and lower surfaces)
thickness = upper_surfacet(:, 2) - lower_surfacet(:, 2);

% Integrate thickness over chord to find cross-sectional area
cross_sectional_area = trapz(upper_surfacet(:,1), thickness);

%volume finding
areat = (cross_sectional_area * chordt.^2 ) / cosd(quarter_chord_sweep);

total_tailplane_volume = (trapz(wing_span1 , areat)) * 2;

figure
hold on
plot(wing_span1 , areat)
title("Area vs Span")
hold off

figure
hold on
plot(wing_span1 , chordt)
title("Area vs Span")
hold off


figure
hold on
plot(upper_surface(:,1) , upper_surface(:,2))
hold off









