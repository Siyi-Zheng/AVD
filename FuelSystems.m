clc
clear

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

disp(["Final CSA for 1m of chord length is" num2str(final_area)])

% Define a range of chord lengths to calculate the CSA
chord_lengths = [0,1.0,1.5,2.0]; % Example chord lengths in meters

% Calculate CSA for each chord length and display
for i = 1:length(chord_lengths)
    chord = chord_lengths(i);
    CSA = final_area * chord^2;
    disp(['For a chord length of ', num2str(chord), ' m, the cross-sectional area is: ', num2str(CSA), ' square meters']);
end

%calculating chord distribution

%leading edge sweep angle
b = 32.45 * tand(26.6);
a = b - 3.47/4;

aa = (a + 3.47)/32.45;

%sweep_angle = atand(aa);
sweep_angle = 29.6;

%chord distribution for first section
outer_chord = 13.89 - (9.75 * tand(sweep_angle));

gradient1 = (outer_chord - 13.89)/9.75;

%for outer section
gradient2 = (3.47 - outer_chord)/(32.45 - 9.75);

span1 = 0:0.01:9.75;
span2 = 9.75:0.01:32.45;

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
area1 = final_area * chord1.^2;
area2 = final_area * chord2.^2;

figure
hold on
plot(span1, area1)
plot(span2, area2)
title("Area Ratio")
hold off

%calculating available wing volume
volume_available1 = 0.5 * (13.89 + outer_chord) * 9.75;
volume_available2 = 0.5 * (outer_chord + 3.47) * (32.45 - 9.75);

total_availablevolume = volume_available1 + volume_available2;

%fuel volume required at 15 degrees density of 804kg/m^3

Wf = 1.760909890249436e+05;

rhof = 804;

volume_hold = Wf / rhof;

volume_holdl = volume_hold * 1000;

%increasing to factor foam area, tank area and extra volume for desnity

volume_req = volume_hold * 1.15;

volume_wings = (volume_hold * 0.60) * 2;

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

fueltank = [upper_surfacea(18:73 , :) ; flipud(lower_surfacea(18:73 , :)) ; upper_surfacea(18 , :)];
%central tank and front tank seperation
%structural member coordinates

% Separate upper and lower surfaces
upper_surface = fueltank(1:56 , :); % Columns 1-46 for upper surface
lower_surface = fueltank(57:112 , :); % Columns 47-92 for lower surface

% Calculate thickness between upper and lower surfaces
thickness = upper_surface(: , 2) - flipud(lower_surface(: , 2));

% Calculate the area by integrating thickness along the chord positions
Tankarea = trapz(upper_surface(: , 1), thickness);

% Display the area
disp(['The total cross-sectional area of the fuel tank is: ', num2str(Tankarea), ' square units']);

% Plot the airfoil
figure
hold on
plot(upper_surfacea(:, 1), upper_surfacea(:, 2), 'b', 'DisplayName', 'Upper Surface');
plot(lower_surfacea(:, 1), lower_surfacea(:, 2), 'r', 'DisplayName', 'Lower Surface');
plot(fueltank(:,1),fueltank(:,2),'k','DisplayName','Fuel Tank')
xlabel('Chord Position (x/c)');
ylabel('Thickness (y/c)');
title('Airfoil Profile');
legend;
ylim([-0.8 0.8]);
hold off

%plotting area distribution
area3 = Tankarea * chord1.^2;
area4 = Tankarea * chord2.^2;

figure
hold on
plot(span1, area3)
plot(span2, area4)
title("Area Ratio")
hold off

%volume calculations
volume1 = (trapz(span1 , area3) * 0.85);
volume2 = trapz(span2 , area4) * 0.85;








