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

%fixed variables
% if there is an error run wing_load_distribution.m first
momentMax = 4362716978.29012/100; % Nm % maximum bending moment
c = (0.7 - 0.12) .* 12.41; % m %wingbox width
b2 = 0.097 * 12.41; % m     web height (wingbox height)  %
E_composite = 120.5 ; %Gpa
E_Aluminium = 73.85; % GPa     %Young's modulus      
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
optimal_rib= Rib_spacing(index)

% spar

a = optimal_rib; % Web panel spacing, m
E = E_Aluminium; % Young's Modulus, GPa
sigma_y = 495; % Yield Stress, MPa
Ks = 7; % Read from graph
V = 3311650; % Shear Load, N
T = 4200000; % Torque Load, Nm
q0 = T / (2 * c * b2 * 1000); % Torque shear flow, N/mm
q2 = V / (2 * b2 * 1000); % Load shear flow, N/mm
q_FS = q2 + q0; % Front spar shear flow, N/mm
q_RS = q2 - q0; % Rear spar shear flow, N/mm
t_FS = (q_FS * 1000 * b2 / (Ks * E * 1e9)) ^ (1/3) * 1000; % Front web thickness, mm
t_RS = (q_RS * 1000 * b2 / (Ks * E * 1e9)) ^ (1/3) * 1000; % Rear web thickness, mm

%compression-shear combination for stringer and skin
tau_0 = q0 / t2; % shear stress, N/mm
b1 = b; % skin panel width, mm
tau_cr = Ks * E_Aluminium * 1e9 * (t2 / b1) ^ 2 * 1e-6; % critical shear buckling stress, MPa
R_c = sigma_0 / sigma_cr; % Compressive stress ratio
R_s = tau_0 / tau_cr; % Shear stress ratio
val = R_s^2 + R_c; % combined stress ratio

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
    disp(loc);
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
Ks = 8.1; % Read from graph
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

for i =1: length(difference_array)
    index = 



figure
plot(span, t_RS)
hold on
plot(span, t_FS)
hold off

