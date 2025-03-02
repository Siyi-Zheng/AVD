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

%fixed variables#
density_comp = 1515;
density_al = 2800;
% if there is an error run wing_load_distribution.m first
momentMax = 4362716978.29012/100; % Nm % maximum bending moment
c = (0.7 - 0.12) .* 12.41; % m %wingbox width
b2 = 0.097 * 12.41; % m     web height (wingbox height)  %
E_composite = 61 ; %Gpa
E_Aluminium = 72.5; % GPa     %Young's modulus      
N = momentMax / (c .* b2); % N/m compressive load per unit length
flangeWebRatio = 0.3; % ratio of flange length to web height

n = 5: 1: 50;   % no of panels
% ts = linspace(1,15, length(n)); % mm      %stringer thickness
% h = linspace(20, 300, length(n)); % mm      %stringer web height
ts= 1:0.1:15;     %stringer thickness
h = 60:2:120;     %stringer web height
b = c ./ n.* 1000; % mm      % 
d = h .* flangeWebRatio; % mm  %flange width of stringer

for i = 1:length(n)
    for j=1:length(ts)
        for k = 1:length(h)
            % TODO: consider buckling of stringers, I think this is why we
            % get unrealistic results with lots of tiny stringers
            As(i,j,k)= ts(j) .* (h(k) + 2 * d(k)); % mm^2
            t2(i,j,k)= (N /3.62/ (E_composite*10^9) * (b(i) / 1000) ^ 2) ^ (1/3) *1000; % mm
            sigma_0(i,j,k) = N / (t2(i,j,k) * 1000); % MPa
            t_e(i,j,k) = t2(i,j,k) + (As(i,j,k) / b(i)); % mm
            effectivePanelCSA(i,j,k) = n(i)* b(i)* t_e(i,j,k); % mm^2
            % catchpole diagram
            catchpoleY(i,j,k) = ts(j)/ t2(i,j,k); % x-value to read
            catchpoleX(i,j,k) = As(i,j,k)/ (b(i)* t2(i,j,k)); % y-value to read
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

% shitty mass optimisation
% F_factor = F_factor ./ As;

sigma_cr = stressRatio .* sigma_0;

%Farrar efficiency (stringer local buckling)
% for g = 1:length(Rib_spacing)


% 3d plot
F_factor(F_factor <= 0.5) = NaN;
for i = 1:length(h)
    slice = squeeze(F_factor(:,:,i));
    % mass = squeeze(effectivePanelCSA(:,:,i));
    xData = n;
    yData = ts;
    color = h(i) .* slice' ./ slice';
    surf(n, ts, slice', color, "EdgeColor", "none");
    % surf(n, ts, mass')
    hold on
end
colormap("turbo")
bar = colorbar;
bar.Label.String = "Stringer Height (mm)";
xlabel("Number of Stringers")
ylabel("Stringer Thickness (mm)")
zlabel("Farrar Efficiency")

% get the optimal values based on the farrar plot

[mxv, idx] = max(F_factor(:));
[id1, id2, id3] = ind2sub(size(F_factor), idx);
n = n(id1);
ts = ts(id2);
h = h(id3);
b = c ./ n.* 1000; % mm


% extract other values at this optimum
As = As(id1, id2, id3);
t2 = t2(id1, id2, id3);   %skin thickness
sigma_0 = sigma_0(id1, id2, id3);
sigma_cr = sigma_cr(id1, id2, id3);

Rib_spacing =  0.5:0.05:7; % m

sigma_f= mxv./sqrt(Rib_spacing./N./(E_Aluminium*10^9))./10^6;
resultant = sigma_cr - sigma_f;
[result, index] = min(abs(resultant));

%optimal rib spacing
optimal_rib= Rib_spacing(index);

% spar

a = optimal_rib; % Web panel spacing, m
E = E_Aluminium; % Young's Modulus, GPa
sigma_y = 495; % Yield Stress, MPa
Ks = 8.5; % Read from graph
V = 3311650; % Shear Load, N
T = 4200000; % Torque Load, Nm
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
F = M^2 * s * b2 * (t_e / 1000) * c / (2 * E_Aluminium * 1e9) / (I ^ 2); % Crush load, N
t_r = 7.3; % design rib thickness, mm
sigma_c = (F / t_r) * (c / 1000); % crush stress, MPa
sigma_b = F / (t_r * c * 1e3); % buckling stress, MPa
t_r2 = F / (sigma_y * c * 1e3); % design rib thickness (for buckling), MPa


%varying rib spacing with bending moments and then get how the rib
%thickness varies with rib spacing.

figure
constant = optimal_rib * momentMax^(1/3);
M_x = max(abs(M_x));
rib_spacing = constant * M_x.^(-1/3); % run the WLD2 code first
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
I = c_span .* (t_e ./ 1000) .^ 3 ./ 12 + c_span .* (t_e ./ 1000) .* b2_span .^ 2 / 4;
crush_load = M_x .^ 2 .* rib_spacing .* b2_span .* (t_e ./ 1000) .* c_span...
    ./ (2 .* E_Aluminium .* 1e9 .* (I .^ 2)) ./ 10;
t_r = (crush_load ./ c_span ./ 3.62 ./ (E_Aluminium .* 1e9) .* b2_span .^ 2) .^ (1/3) .* 1000;
rib_thickness = max(minThickness, interp1(span, t_r, rib_list));
% % calculation
% chord = chord * 0.58;
% I = chord .* (t_e ./ 1000) .^ 3 ./ 12 + chord .* (t_e / 1000) * b2 .^ 2 / 4;
% crush_load = M_x .^ 2 .* rib_spacing .* b2 .* (t_e ./ 1000) .* chord...
%     ./ (2 .* E_Aluminium .* 1e9 .* (I .^ 2)) ./ 10;
% t_r = (crush_load ./ chord ./ 3.62 ./ (E_Aluminium .* 1e9) .* b2 .^ 2) .^ (1/3) .* 1000;
% rib_thickness = max(minThickness, interp1(span, t_r, rib_list));
scatter(rib_list, rib_thickness, "kx");
area_list = 0.0939 * c_span .^2;
rib_area_list = [];
for i = 1:length(rib_list)
    rib_loc = rib_list(i);
    rib_area_list = [rib_area_list max(area_list(span > rib_loc))];
end
main_rib_mass = sum(rib_area_list .* rib_thickness / 1000) * density_al;
xlabel("Semi-span (m)")
ylabel("Rib Thickness (mm)")
ylim([0 max(rib_thickness) + 1])
grid on

% plotting how the web thickness varies with the span 
% IF RIB SPACING CHANGES THEN KS CHANGES
% DOUBLE CHECK KS VALUES
V_x = max(abs(shear_x));
V_x(721) = [];
T = max(abs(M_torque));
T(721) = [];
a = rib_thickness; % Web panel spacing, m
E = E_Aluminium; % Young's Modulus, GPa
sigma_y = 495; % Yield Stress, MPa
Ks = 8.5; % Read from graph
V = V_x; % Shear Load, N
q0 = T ./ (2 .* c_span .* b2_span .* 1000); % Torque shear flow, N/mm
q2 = V ./ (2 .* b2_span .* 1000); % Load shear flow, N/mm
q_FS = abs(q2 + q0); % Front spar shear flow, N/mm
q_RS = abs(q2 - q0); % Rear spar shear flow, N/mm
t_FS = (q_FS .* 1000 .* b2_span ./ (Ks .* E_Aluminium .* 1e9)) .^ (1/3) .* 1000; % Front web thickness, mm
t_RS = (q_RS .* 1000 .* b2_span ./ (Ks .* E_Aluminium .* 1e9)) .^ (1/3) .* 1000; % Rear web thickness, mm

discrete_thickness = 1;
difference_array = abs(t_RS - t_RS(1));
number = max(difference_array) - min(difference_array);
difference_array_f = abs(t_FS - t_FS(1));
number_f = max(difference_array_f) - min(difference_array_f);

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
        t_RS_dis(index(i):index(i+1)-1) = t_RS(index(i));
    else
        t_RS_dis(index(i): length(t_RS)) = t_RS(index(i));
    end
end

index_f = [];  % Initialize index array
for i = 1:number
    idx_f = find(difference_array_f < i & difference_array_f > i - 1, 1);  % Find first occurrence
    if ~isempty(idx_f)
        index_f(end+1) = idx_f;  % Append found index
    end
end

index_f(1) =1;

for i = 1:length(index_f)
    if i < length(index_f)
        t_FS_dis(index_f(i):index_f(i+1)-1) = t_FS(index_f(i));
    else
        t_FS_dis(index_f(i):length(t_FS)) = t_FS(index_f(i));
    end
end

t_FS_dis(t_FS_dis<2) = 2;
t_RS_dis(t_RS_dis<2) = 2;

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
b_span = c_span ./ n.* 1000;% skin panel width, mm
N_span = moment_span ./ (c_span .* b2_span);
t2_span= (N_span ./3.62 ./ (E_composite.*10^9) .* (b_span./ 1000) .^ 2) .^ (1/3) .*1000; % mm

difference_array_new = abs(t2_span - t2_span(1));
number_new = max(difference_array_new) - min(difference_array_new);

index = [];  % Initialize index array
for i = 1:number_new
    idx = find(difference_array_new < i & difference_array_new > i - 1, 1);  % Find first occurrence
    if ~isempty(idx)
        index(end+1) = idx;  % Append found index
    end
end


t2_span(t2_span<1) = 1;

%compression-shear combination for stringer and skin
tau_0 = q0 ./ t2_span; % shear stress, N/mm
tau_cr = Ks .* E_Aluminium .* 1e9 .* (t2_span ./ b_span) .^ 2 .* 1e-6; % critical shear buckling stress, MPa
sigma_0_span= N_span ./ (t2_span .* 1000); % MPa
Stress_Ratio_optimal = interp2(Y_grid, X_grid , stressFactors, ts/t2, As/(b*t2));
sigma_cr_span = Stress_Ratio_optimal .* sigma_0_span;
R_c = sigma_0_span ./ sigma_cr_span; % Compressive stress ratio
R_s = tau_0 ./ tau_cr; % Shear stress ratio
val = R_s.^2 + R_c; % combined stress ratio

figure
plot(span, val)

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
            t_s = sqrt(tau_cr .* 1e6 .* b^2 ./ (buckling_coeff .* E_composite .* 1e9)) .* 1000;
        end
        skin_mass = skin_mass + b .* a .* t_s./1000 .* density_comp;
        rib_mass = rib_mass + density_al * area_list(rib) * rib_thickness / 1000; % kg
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
