% skin
WLD2;
table = readtable("catchpole_data.csv");
catchpoledata = table(1:8, 1:10);
farrardata= table(9:16, 1:10);

% % max. bending moment - no stringers

% t2 = (N ./ (3.62 * E ./ c .^ 2)) .^ (1/3) *1000; % mm
% sigma = N ./ (t2 * 1000); % MPa Critical buckling stress
% panelCSA = c .* t2 * 1000; % mm^2


%%%%%%%%%Design panel size and stringer size
% max bending moment - with stringers

%material properties
density_T861 = 2765; % all kg/m^3
density_T36 = 2765;
density_H34 = 2655;
sigma_y_T861 = 431; % Yield Stress, MPa
E_T861 = 73.85; % all GPa
E_T36 = 73.85;
E_H34 = 70.3;

% if there is an error run wing_load_distribution.m first
momentMax = 4362716978.29012/100; % Nm % maximum bending moment
c = (0.7 - 0.12) .* 12.4085; % m %wingbox width
b2 = 0.097 * 12.4085; % m     web height (wingbox height)  %
N = momentMax / (c .* b2); % N/m compressive load per unit length
flangeWebRatio = 0.3; % ratio of flange length to web height

n = 5: 1: 48;   % no of panels
% ts = linspace(1,15, length(n)); % mm      %stringer thickness
% h = linspace(20, 300, length(n)); % mm      %stringer web height
ts= 1:0.1:15;     %stringer thickness
h = 60:2:120;     %stringer web height
b = c ./ n.* 1000; % mm      % Panel width
d = h .* flangeWebRatio; % mm  %flange width of stringer
rib_spacing = 0.5; %vary this after iterations
for i = 1:length(n)
    for j=1:length(ts)
        for k = 1:length(h)
            As(i,j,k)= ts(j) .* (h(k) + 2 * d(k)); % mm^2
            t2(i,j,k)= (N /3.62/ (E_T861*10^9) * (b(i) / 1000) ^ 2) ^ (1/3) *1000; % mm
            sigma_0(i,j,k) = N / (t2(i,j,k) * 1000); % MPa
            t_e(i,j,k) = t2(i,j,k) + (As(i,j,k) / b(i)); % mm
            effectivePanelCSA(i,j,k) = n(i)* b(i)* t_e(i,j,k); % mm^2
            % catchpole diagram
            catchpoleY(i,j,k) = ts(j)/ t2(i,j,k); % x-value to read
            catchpoleX(i,j,k) = As(i,j,k)/ (b(i)* t2(i,j,k)); % y-value to read
             %%optimise using mass
            volume(i,j,k) = n(i) * b(i) * t_e(i,j,k) * rib_spacing *1000;
            mass(i,j,k) = ((density_T861/(1000)^3) * volume(i,j,k));
        end
    end
end

catchpoleY(catchpoleY > 1.4 | catchpoleY < 0.4) = 0;
catchpoleX(catchpoleX>1.2) = 0;
Y_grid = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4]; % 1x10 (Y-coordinates)
X_grid = [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2];  % 8x1 (X-coordinates), NaN ignored
catchpoledata_matrix = table2array(catchpoledata);
farrardata_matrix = table2array(farrardata);
stressFactors = catchpoledata_matrix(2:end,2:end);
f_testing = farrardata_matrix(2:end,2:end);

for i = 1:length(n)
    for j=1:length(ts)
        for k = 1:length(h)
            if catchpoleY(i,j,k) ==0 || catchpoleX(i,j,k) == 0
                continue
            end
            stressRatio(i,j,k) = interp2(Y_grid,X_grid, stressFactors, catchpoleX(i,j,k),catchpoleY(i,j,k)); 
            F_factor(i,j,k) = interp2(Y_grid,X_grid, f_testing, catchpoleX(i,j,k),catchpoleY(i,j,k));
        end
    end
end

% F_factor = F_factor ./ As;

sigma_cr = stressRatio .* sigma_0;

% get the optimal values based on the farrar plot
[if1, if2, if3] = ind2sub(size(F_factor), find(F_factor>0.7));
% Extract values from mass at the found indices
selectedValues = mass(sub2ind(size(mass), if1, if2, if3));

% Find the minimum value and its index
[minValue, minIdx] = min(selectedValues);

% % 3d plot
figure
F_factor(F_factor <= 0.7) = NaN;
for i = 1:length(h)
    slice = squeeze(mass(:,:,i));
    % mass = squeeze(effectivePanelCSA(:,:,i));
    xData = n;
    yData = ts;
    color = h(i) .* slice' ./ slice';
    surf(n, ts, slice', color, "EdgeColor", "none");
    % surf(n, ts, mass')
    hold on
end
X= xData(if1);
Y= yData(if2);
Z = selectedValues;
scatter3(X, Y, Z, "x")


colormap("turbo")
bar = colorbar;
bar.Label.String = "Stringer Height (mm)";
xlabel("Number of Stringers")
ylabel("Stringer Thickness (mm)")
zlabel("Mass")

% Get the corresponding row, column, and depth indices
id1 = if1(minIdx);
id2 = if2(minIdx);
id3 = if3(minIdx);
% new_mass= mass([if1, if2, if3])
% [mxv, idx] = min(mass(:));
% [id1, id2, id3] = ind2sub(size(F_factor), idx);
n = n(id1)
ts = ts(id2)
h = h(id3)
b = c ./ n.* 1000; % mm
mxv = F_factor(id1, id2, id3)


% extract other values at this optimum
As = As(id1, id2, id3);
t2 = t2(id1, id2, id3)  ; %skin thickness
sigma_0 = sigma_0(id1, id2, id3);
sigma_cr = sigma_cr(id1, id2, id3);

RATIO1= As/(b*t2)
RATIO2 = ts / t2

Rib_spacing =  0.5:0.05:2; % m

sigma_f= mxv./sqrt(Rib_spacing./N./(E_T861*10^9))./10^6;
resultant = sigma_cr - sigma_f;
[result, index] = min(abs(resultant));

%optimal rib spacing
optimal_rib_idk= Rib_spacing(index)

optimal_rib = 0.5;

% spar

a = optimal_rib; % Web panel spacing, m
E = E_H34; % Young's Modulus, GPa
Ks = 8.5; % Read from graph
V = 3311650; % Shear Load, N
T = 2.5e6; % Torque Load, Nm
q0 = T / (2 * c * b2 * 1000); % Torque shear flow, N/mm
q2 = V / (2 * b2 * 1000); % Load shear flow, N/mm
q_FS = q2 + q0; % Front spar shear flow, N/mm
q_RS = q2 - q0; % Rear spar shear flow, N/mm
t_FS = (q_FS * 1000 * b2 / (Ks * E * 1e9)) ^ (1/3) * 1000; % Front web thickness, mm
t_RS = (q_RS * 1000 * b2 / (Ks * E * 1e9)) ^ (1/3) * 1000; % Rear web thickness, mm

% ribs
M = momentMax;
t_e = ts + (As / b);
s = 2; % rib spacing, m
I = (c * (t_e / 1000) ^ 3 / 12 + c * (t_e / 1000) * (b2 / 2) ^ 2); % 2nd moment of area of panel, m^4
F = M^2 * s * b2 * (t_e / 1000) * c / (2 * E_T861 * 1e9) / (I ^ 2); % Crush load, N
t_r = 7.3; % design rib thickness, mm
sigma_c = (F / t_r) * (c / 1000); % crush stress, MPa
sigma_b = F / (t_r * c * 1e3); % buckling stress, MPa
t_r2 = F / (sigma_y_T861 * c * 1e3); % design rib thickness (for buckling), MPa

%varying rib spacing with bending moments and then get how the rib
%thickness varies with rib spacing.

figure
constant = optimal_rib * momentMax^(1/3);
M_x = max(abs(M_x));
rib_spacing = constant * M_x.^(-1/3) .* (0.5/0.495); % run the WLD2 code first
rib_list = [min(span)];
while i < length(span)
    loc = rib_list(end);
    next_rib = loc + rib_spacing(i);
    i = find(span > next_rib, 1);
    rib_list = [rib_list next_rib];
end
rib_list = rib_list(1:end-1);

% rib thickness calculation
% remove duplicates
minThickness = 2; % mm
span = unique(span);
M_x(721) = [];
chord = flip(unique(chord));
rib_spacing(721) = [];
b2_span = 0.097 * chord;
c_span = chord * 0.58;
I = c_span .* (t_e ./ 1000) .^ 3 / 12 + c_span .* (t_e ./ 1000) .* b2_span .^ 2 / 4;

crush_load = M_x .^ 2 .* rib_spacing .* b2_span .* (t_e ./ 1000) .* c_span...
    ./ (2 .* E_T861 .* 1e9 .* (I .^ 2)) ./ 10;
t_r = (crush_load ./ c_span ./ 3.62 ./ (E_T861 .* 1e9) .* b2_span .^ 2) .^ (1/3) .* 1000;
rib_thickness = max(minThickness, interp1(span, t_r, rib_list));
% % calculation
% chord = chord * 0.58;
% I = chord .* (t_e ./ 1000) .^ 3 ./ 12 + chord .* (t_e / 1000) * b2 .^ 2 / 4;
% crush_load = M_x .^ 2 .* rib_spacing .* b2 .* (t_e ./ 1000) .* chord...
%     ./ (2 .* E_T861 .* 1e9 .* (I .^ 2)) ./ 10;
% t_r = (crush_load ./ chord ./ 3.62 ./ (E_T861 .* 1e9) .* b2 .^ 2) .^ (1/3) .* 1000;
% rib_thickness = max(minThickness, interp1(span, t_r, rib_list));
scatter(rib_list, rib_thickness, 70, "kx");
area_list = 0.0939 * c_span .^2;
rib_area_list = [];
for i = 1:length(rib_list)
    rib_loc = rib_list(i);
    rib_area_list = [rib_area_list max(area_list(span > rib_loc))];
end
main_rib_mass = sum(rib_area_list .* rib_thickness / 1000) * density_T861;
xlabel("Semi-span (m)")
ylabel("Rib Thickness (mm)")
ylim([0 max(rib_thickness) + 1])
grid on


%mass optimisation to get optimal rib spacing
rib_spacing_test = 0.2:0.05:1.5; % m %vary this after iterations

number_ribs = (32.5 ./rib_spacing_test) +1;
t_e_farrar = (1/mxv) .* sqrt(N .* rib_spacing_test./ E_T861); %mm

volume= n .* b .* (t_e_farrar) .* rib_spacing_test.*1000; %mm^3
mass_skin_stringer = (((density_T861.* 1e-9) .* volume))  ./(span(end) * c);

mass_ribs = number_ribs .* ((rib_area_list(1).* rib_thickness(1) ./ 1000) .* density_T861) ./(span(end) * c);

total_mass = (mass_skin_stringer + mass_ribs);

figure 
plot(rib_spacing_test, total_mass, Linewidth = 1.5)
hold on
plot(rib_spacing_test, mass_skin_stringer, Linewidth = 1.5)
plot(rib_spacing_test, mass_ribs, Linewidth = 1.5)
xlabel("Rib spacing")
ylabel("Mass per unit area (Kg/m^2)")
grid on
legend('Total', 'Skin-stringer', 'Ribs')
hold off

[value, in ]= min(total_mass);
final_optimal_rib_spacing = rib_spacing_test(in);

% plotting how the web thickness varies with the span 
% IF RIB SPACING CHANGES THEN KS CHANGES
% DOUBLE CHECK KS VALUES
V_x = max(abs(shear_x));
V_x(721) = [];
T = max(abs(M_torque));
T(721) = [];
a = rib_thickness; % Web panel spacing, m
E = E_T861; % Young's Modulus, GPa
Ks = 8.1; % Read from graph
V = V_x; % Shear Load, N
q0 = T ./ (2 .* c_span .* b2_span .* 1000); % Torque shear flow, N/mm
q2 = V ./ (2 .* b2_span .* 1000); % Load shear flow, N/mm
q_FS = abs(q2 + q0); % Front spar shear flow, N/mm
q_RS = abs(q2 - q0); % Rear spar shear flow, N/mm
t_FS = (q_FS .* 1000 .* b2_span ./ (Ks .* E_H34 .* 1e9)) .^ (1/3) .* 1000; % Front web thickness, mm
t_RS = (q_RS .* 1000 .* b2_span ./ (Ks .* E_H34 .* 1e9)) .^ (1/3) .* 1000; % Rear web thickness, mm

t_FS_dis = discretisation_1mm(t_FS, 2);
t_RS_dis = discretisation_1mm(t_RS, 2);

figure
plot(span, t_RS, LineStyle="--", LineWidth= 1.5)
hold on
plot(span, t_FS, LineStyle="--", LineWidth= 1.5)
plot(span, t_FS_dis, LineWidth=1.5)
plot(span, t_RS_dis, LineWidth=1.5)
hold off
grid on
xlabel("Semi span")
ylabel("Spar web thickness (mm)")
legend("Theoretical rear spar thickness", "Theoretical front spar thickness", "Discretised front spar thickness", "Discretised rear spar thickness")

%%%%%%%%%shear flow in skin
moment_span = M_x;
% b_span = c_span ./ n.* 1000;% skin panel width, mm
b_span = b;
N_span = moment_span ./ (c_span .* b2_span);
t2_span= (N_span ./3.62 ./ (E_T861.*10^9) .* (b_span./ 1000) .^ 2) .^ (1/3) .*1000; % mm
t2_span(t2_span<1) =  1;
t2_span_disc = discretisation_1mm(t2_span, 1);

%compression-shear combination for stringer and skin
tau_0 = q0 ./ t2_span; % shear stress, N/mm
tau_cr = Ks .* E_T861 .* 1e9 .* (t2_span./ b_span) .^ 2 .* 1e-6; % critical shear buckling stress, MPa
sigma_0_span= N_span ./ (t2_span .* 1000); % MPa
Stress_Ratio_optimal = interp2(Y_grid, X_grid , stressFactors, ts./t2_span, As./(b.*t2_span));
sigma_cr_span = Stress_Ratio_optimal .* sigma_0_span;
R_c = sigma_0_span ./ sigma_cr_span; % Compressive stress ratio
R_s = tau_0 ./ tau_cr; % Shear stress ratio
val = R_s.^2 + R_c; % combined stress ratio

disp(max(val))

figure
plot(span, val, LineWidth=1.5)
grid on
xlabel("Semi-span (m)")
ylabel("R_s ^2 + R_c")

% Get the final mass of all components
t_e_span = t2_span_disc + As ./ b;
[unique_t_e, stations] = unique(t_e_span);
stations_r = flip(stations);
span_ss= span(stations_r);
l_y_ss= (0.7-0.12).*chord(stations_r);


for jn=1:length(l_y_ss)
    if jn == length(l_y_ss)
        area_ss(jn)= 0.5*(l_y_ss(jn) + (0.7-0.12)*chord(end)) *(span(end)-span_ss(jn));
    else
        area_ss(jn) = 0.5* (l_y_ss(jn) + l_y_ss(jn+1)) * (span_ss(jn+1)-span(jn));
    end
end

volume_ss = sum(area_ss.* (unique_t_e./1000)); %m^3
skin_stringer_fmass = ((density_T861) .* volume_ss)*2;

% n_span = floor(c_span(stations_r) ./ (b) .*1000)
% dx= [diff(span_ss), span(end)-span_ss(end)];
% y_ss = (1-0.41).* chord(stations_r)
% length_y_ss = [diff(y_ss), y_ss(end)- (1-0.41)*chord(end)]
% length_ss = sqrt(dx.^2 + length_y_ss.^2)
% volume_ss= n_span.* b .* unique_t_e .*1000 .* length_ss.*1000 %mm^3
% masse= (density_T861.* 1e-9) .* volume_ss
% skin_stringer_fmass = (sum(((density_T861.* 1e-9) .* volume_ss)))*2;

ribs_fmass = main_rib_mass;

[unique_fs, ln]= (unique(t_FS_dis));
length_ln = flip(span(ln));
kh = (1-0.12)*chord(ln);
b2_span_unique = b2_span(ln);
length_x= [diff(length_ln), span(end)-length_ln(end)];
length_y = [diff(kh), kh(end)-(1-0.12)*chord(end)];
length_s = sqrt(length_x.^2 + length_y.^2);

volume_FS = unique_fs .* b2_span_unique.*1000 .* length_s.*1000; %mm^3
spar_fs_mass = sum((volume_FS).*(density_H34.* 1e-9));

[unique_rs, lnr]= (unique(t_RS_dis));
length_lnr = flip(span(lnr));
b2_span_uniquer = b2_span(lnr);
khr = (1-0.7)*chord(lnr);
lengthr_x = [diff(length_lnr), span(end)-length_lnr(end)];
lengthr_y = [diff(khr), khr(end)-(1-0.7)* chord(end)];
lengthr = sqrt(lengthr_x.^2 + lengthr_y.^2);
volume_RS = unique_rs .* b2_span_uniquer .*1000 .* lengthr.*1000;
spar_rs_mass = sum((volume_RS).*(density_H34.* 1e-9));

spar_mass = spar_fs_mass + spar_rs_mass;

mass_fuck = ribs_fmass + skin_stringer_fmass + spar_mass
disp(mass_fuck)
%% 

span_old_var = span; % so it doesn't get overwritten

% get cylindrical buckling data ready to interpolate
bucklingdata = readtable("buckling_dcell.csv");
Y_grid = [0, 2, 4, 6, 10, 16, 50]; % (Y-coordinates)
X_grid = [0, 0.2, 0.667, 1, 1.5, 2, 3, 10];  % (X-coordinates)
bucklingdata_matrix = table2array(bucklingdata);
buckling_coeffs = bucklingdata_matrix(2:end,2:end);

% d-cell design
perimeter = 0.1739 * c_span;
area = 0.243 / 100 * c_span .^ 2;
min_pseudoribs = 8;
max_pseudoribs = 45;
root_radius = 0.3; % m, this can be anything lol
radius = root_radius .* c_span / c_span(1);
buckling_coeff = [];
total_mass = [];
rib_thickness = 2; % mm, idk what a sensible value for this is
tau_cr = 145; % tresca shear yielding stress
for num_pseudoribs = min_pseudoribs:max_pseudoribs
    span = span_old_var;
    a = (max(span) - min(span)) / (num_pseudoribs + 1);
    pseudorib_locs = linspace(min(span), max(span), num_pseudoribs);
    % remove outer rib in calculations of d-cell perimeter
    b_list = [];
    radius_list = [];
    span_list = [];
    area_list = [];
    for i = 1:num_pseudoribs
        p = max(perimeter(span >= pseudorib_locs(i)));
        r = max(radius(span >= pseudorib_locs(i)));
        s = max(span(span >= pseudorib_locs(i)));
        ar = max(area(span >= pseudorib_locs(i)));
        b_list = [b_list p];
        radius_list = [radius_list r];
        span_list = [span_list s];
        area_list = [area_list ar];
    end
    aspect_ratio = a./b_list;
    skin_mass = 0; % start running total
    rib_mass = 0;
    for rib = 1:num_pseudoribs
        t_s = 10; % iterative process
        old_ts = 0;
        radius_at_rib = radius_list(rib);
        b = b_list(rib);
        span = span_list(rib);
        bucklingX = aspect_ratio(rib); % x-coordinate of graph to read from
        step = 0; % stops solution diverging
        while abs(old_ts - ts) > 1e-3 && step < 10
            step = step + 1;
            old_ts = t_s;
            bucklingY = min(a, b) / sqrt(radius_at_rib * t_s / 1000); % convert to m
            % actually do the interpolation
            buckling_coeff = interp2(Y_grid, X_grid, buckling_coeffs,...
                bucklingY, bucklingX);
            % shear stress calculation
            t_s = sqrt(tau_cr .* 1e6 .* b^2 ./ (buckling_coeff .* E_T861 .* 1e9)) .* 1000;
        end
        skin_mass = skin_mass + b .* a .* t_s./1000 .* density_T861;
        rib_mass = rib_mass + density_T861 * area_list(rib) * rib_thickness / 1000; % kg
    end
    total_skin_mass(num_pseudoribs) = skin_mass;
    total_rib_mass(num_pseudoribs) = rib_mass;
    total_mass(num_pseudoribs) = total_rib_mass(num_pseudoribs) + total_skin_mass(num_pseudoribs);
end

figure
ax = gca();
plot(8:max_pseudoribs, total_mass(8:end), "k-")
hold on
plot(8:max_pseudoribs, total_skin_mass(8:end), "r-")
ylabel("Skin/Total mass (kg)")
yyaxis right
plot(8:max_pseudoribs, total_rib_mass(8:end), "b-")
grid on
xlim([8, max_pseudoribs])
% Set the color of each axis to black
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = [0 0 0];
xline(25, "k--")
% labels
xlabel("Number of pseudoribs")
ylabel("Total pseudorib mass (kg)")
legend(["Total", "Skin", "Pseudoribs"], Location="northwest")



function disc_val = discretisation_1mm( t_RS, min_val)

difference_array = abs(t_RS - t_RS(1));
number = max(difference_array) - min(difference_array);

index = [];  % Initialize index array
for i = 1:number
    idx = find(difference_array < i & difference_array > i - 1, 1);  % Find first occurrence
    if ~isempty(idx)
        index(end+1) = idx;  % Append found index
    end
end

index(1) =1;

for i = 1:length(index)
    if i < length(index)
        disc_val(index(i):index(i+1)-1) = t_RS(index(i));
    else
        disc_val(index(i): length(t_RS)) = t_RS(index(i));
    end
end

disc_val(disc_val< min_val) = min_val;

end
