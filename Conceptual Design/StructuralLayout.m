clear
clc

% parameters
area = 481.77; % m^2
wing_span = 65; % m
taper_ratio = 0.25;
quarter_chord_sweep = 26.6; % degrees
trailing_edge_kink = 0.3; % fraction of wing span
area_ratio = 0.09161; % cross-sectional area as a proportion of c^2

% simple design
total_chord = 1 + 1 / taper_ratio;
tip_chord = 2 * area / (wing_span * total_chord);
root_chord = tip_chord / taper_ratio;
qcy = -wing_span / 2 * tand(quarter_chord_sweep);
le_tip = qcy + tip_chord / 4;
te_tip = qcy - 3 * tip_chord / 4;
le_root = root_chord / 4;
te_root = -3 * root_chord / 4;



% wing with trailing edge kink
kink_span = wing_span / 2 * trailing_edge_kink;
% recalculate the wing iteratively
A_tot = 0; % initialize the total area to a value larger than the target area
tip_chord2 = tip_chord;
root_chord2 = root_chord;

% iterative loop
while A_tot < area/2
    A_le = wing_span / 16 * (tip_chord2 + root_chord2);
    A_inner = kink_span / 2 * (3/2 * root_chord2 - kink_span * tand(quarter_chord_sweep));
    A_outer = (wing_span / 2 - kink_span)/2 * (3/4 * root_chord2 - kink_span * tand(quarter_chord_sweep) + 3/4 * tip_chord2);
    A_tot = A_le + A_inner + A_outer;
    root_chord2 = root_chord2 + 1e-3;
    tip_chord2 = root_chord2 * taper_ratio;
end

qcy = -wing_span / 2 * tand(quarter_chord_sweep);
le_tip = qcy + tip_chord2 / 4;
te_tip = qcy - 3 * tip_chord2 / 4;
le_root = root_chord2 / 4;
te_root = -3 * root_chord2 / 4;

le_tip = le_root - ( 9.75/cosd(30.16) + 23.1) * sind(30.16);

te_angle = acosd(18824.625/19978.13);

length = ((wing_span/2) - 9.75 - 18.824625);
length1 = 18.82625 + length;

te_tip = te_root - sind(te_angle) * length1;


grad_le = (le_tip-le_root)/((wing_span)/2);
grad_te = (te_tip-te_root)/(((wing_span)/2)-kink_span);

x_te = 0:0.01:(wing_span/2);
x_le = 0:0.01:((9.75/cosd(30.16) + 23.1) * cosd(30.16));

y_le = x_le * grad_le + le_root;
y_te = (x_te(1:find(x_te == 9.75)) * 0) + te_root;
c = te_root - grad_te * kink_span;
y_te = ([y_te , (x_te(find(x_te == 9.76):end)) * grad_te + c]);


index_box = 256;
y_lebox = y_le;
y_tebox = y_te;
y_lebox(1:index_box) = y_lebox(index_box);
y_tebox(1:index_box) = y_tebox(index_box);

%%
%adding le raked wingtip
x_vals = [x_le(end) , 31  , 32  , x_te(end) ];
y_vals = [ y_lebox(end) , -13.5  , -16  , y_tebox(end)];

% Define a finer set of x-values for smooth evaluation
x_fine = x_le(end):0.01:x_te(end); 

end_gradients = [grad_le , -10];

% Fit a cubic spline to the control points
spline_fit = spline(x_vals, y_vals);

% Evaluate the spline at the finer x-values
y_fine = ppval(spline_fit, x_fine);


%%
%adding central box
x_box = -2.55:0.01:2.55;
y_leboxa = (x_box * 0) + 2.55 * grad_le + le_root;
y_teboxa = (x_box * 0) + te_root;


%%
%adding spars
%le_spar @ 12% chord
chord = y_lebox - y_tebox(1:2973);
front_spar = y_lebox - 0.12 * chord;
rear_spar = y_lebox - 0.7 * chord;


%extending spars
%rear spar extension
grad_rspar = (rear_spar(end) - rear_spar(976)) / (x_le(end) - 9.75);
grad_fspar = (front_spar(end) - front_spar(976)) / (x_le(end) - 9.75);

cr = rear_spar(976) - grad_rspar * 9.75;
cf = front_spar(976) - grad_fspar * 9.75;

x_extension = x_le(end)+0.01:0.01:x_te(end);

front_spar = [front_spar , (grad_fspar * x_extension + cf)];
rear_spar = [rear_spar , (grad_rspar * x_extension + cr)];





%%
%plotting heavy ribs
%ribs to reinforce wing box
rib_6 = [0.5 , 0.5 ; y_lebox(51) , y_tebox(51)];

%rib one at 2.55
rib_1 = [2.55 , 2.55 ; y_le(256) , y_te(256)];

%rib 2 at undercarriage mounting to transmit load
rib_2 = [5.6 , 5.6 ; y_le(561) , y_te(561)];

%rib 3 at 30% span
find = 0.3 * wing_span/2;
%index = find(x_le == find);
rib_3 = [find , find ; y_le(976) , y_te(976)];



%rib 4 at 65% span
find = round(0.65 * wing_span/2 , 2);
index = 2114;
rib_4 = [find , find ; y_le(index) , y_te(index)];

%rib 5 at end of spars
rib_5 = [x_le(end) , x_le(end) ; y_le(end) , y_te(2973)];



%%
%Flaps
%front line of flaps is 25% chord, and flaps range from 10% span to 65%
span_64 = 0.64 * wing_span/2;
span_10 = 0.1 * wing_span/2;

index1 = round(span_10 * 100) + 2;
index2 = round(span_64 * 100) - 1;

flap_le = y_le(index1:index2) - 0.75 .* chord(index1:index2);
flap_te = y_te(index1:index2);

FLAP = [ x_le(index1:index2) , fliplr(x_le(index1:index2)) , x_le(index1) ; flap_le , fliplr(flap_te) , flap_le(1) ];

%split flap to inner and outer
inner_flap = [FLAP(:,1:635) , FLAP(:, 2873:end)];
outer_flap = FLAP(:,665:2843);

%%
%adding aeloerons
span_67 = 0.67 * wing_span/2;
aileron_end = 2960;

index3 = round(span_67 * 100);

aileron_le = y_le(index3:aileron_end) - 0.75 .* chord(index3:aileron_end);
aileron_te = y_te(index3:aileron_end);

AILERON = [x_le(index3:aileron_end), fliplr(x_le(index3:aileron_end)) , x_le(index3) ; aileron_le , fliplr(aileron_te) , aileron_le(1)];

%%
%adding slats

slat_le = y_le(index1:end);
slat_te = y_le(index1:end) - 0.115 .* chord(index1:end);

SLAT = [x_le(index1:end) , fliplr(x_le(index1:end)) , x_le(index1) ; slat_le , fliplr(slat_te) , slat_le(1)];

%split flap to inner and outer
inner_slat = [SLAT(:,1:635) , SLAT(:, 4660:end)];

centre_slat = [SLAT(:,665:1773) , SLAT(:, 3522:4630 )];

outer_slat = [SLAT(: , 1803:2632) , SLAT(: , 2663:3492)];


%%

%fuel tanks
tank_le = front_spar(261:2973) - 0.1;
tank_te = rear_spar(261:2973) + 0.1;

TANK = [x_le(261:end) , fliplr(x_le(261:end)) , x_le(261) ; tank_le , fliplr(tank_te) , tank_le(1)];

index_inc = index1 - 261;

%split flap to inner and outer
inner_tank = [TANK(:,1:701) , TANK(:, 4726:end)];

centre_tank = [TANK(:,731:1839) , TANK(:, 3588:4696 )];

outer_tank = [TANK(: , 1869:2698) , TANK(: , 2729:3558)];

%main fuselage tank
length = 8.72;

front = -3.3 + length;

fuselage_tank = [-2.6 , 2.6 , 2.6 , -2.6 ; -3.3 , -3.3 , front , front] ;

%%
%undercarriage
L5 = 2.992;
L6 = 1.451;
theta = 14; %degrees
L1 = 3.388815;
L2 = 0.695146;
L3 = 1.467443;
L4 = 2.954526;


p1 = [L6; (front_spar(1) - L5)];
p2 = [p1(1) + cosd(14) * L3; p1(2) - sind(14) * L3];
temp = (L4 - L2) / 2;
p3 = [p2(1) - sind(14) * temp; p2(2) - cosd(14) * temp];
p4 = [p3(1) + cosd(14) * L1 ; p3(2) - sind(14) * L1];
p5 = [p4(1) - sind(14) * L2; p4(2) - cosd(14) * L2];
p6 = [p3(1) - sind(14) * L2; p3(2) - cosd(14) * L2];
p7 = [p6(1) - sind(14) * temp; p6(2) - cosd(14) * temp];
p8 = [p7(1) - cosd(14) * L3; p7(2) + sind(14) * L3];

p1(2) = p1(2) - 0.95;
p2(2) = p2(2) - 0.95;
p3(2) = p3(2) - 0.95;
p4(2) = p4(2) - 0.95;
p5(2) = p5(2) - 0.95;
p6(2) = p6(2) - 0.95;
p7(2) = p7(2) - 0.95;
p8(2) = p8(2) - 0.95;

%creating box for undercarriage
points = [p1, p2, p3, p4, p5, p6, p7, p8, p1]; % Include p1 again to close the loop

%initialising
x_values = [];
y_values = [];

for i=2:9

    grad = (points(2,i) - points(2,i-1))/(points(1,i) - points(1,i-1));
    c = points(2,i) - grad * points(1,i);

    if points(1, i-1) < points(1, i)
        x = points(1, i-1):0.01:points(1, i);
    else
        x = points(1, i-1):-0.01:points(1, i);
    end
    y = grad .* x + c;

    x_values = [x_values , x];
    y_values = [y_values , y];

end




% Filter for x > 2.55
x_filtered = x_values(x_values > 2.60);
y_filtered = y_values(x_values > 2.60);

vals = [fliplr(x_filtered) ; fliplr(y_filtered)];

%add to tank values
inner_tank = [inner_tank , vals];

%%
%finally plotting fuselage lines

fuselage = [3.17 , 3.17 ; 50 , -50];


%% PLOTTING

figure
figure('Position', [80, 80, 1200, 400]); 

hold on

% Plotting wing
plot(x_le(256:end), y_le(256:end), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Wing Planform'); % Included in legend
plot(x_te(256:end), y_te(256:end), 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off'); % Excluded from legend
plot(-x_le(256:end), y_le(256:end), 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(-x_te(256:end), y_te(256:end), 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(x_fine, y_fine , 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off')
plot(-x_fine, y_fine , 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off')

% Plotting wing box
plot(x_box, y_leboxa, 'Color', [0.9290 0.6940 0.1250], 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'WingBox'); % Included in legend
plot(x_box, y_teboxa, 'Color', [0.9290 0.6940 0.1250], 'LineStyle', '--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot([-2.55, -2.55], [y_leboxa(1), y_teboxa(1)], 'Color', [0.9290 0.6940 0.1250], 'LineStyle', '--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot([2.55, 2.55], [y_leboxa(end), y_teboxa(end)], 'Color', [0.9290 0.6940 0.1250], 'LineStyle', '--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Plotting spars
plot(x_te(1:3101), front_spar(1:3101), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Spar'); % Included in legend
plot(x_te(1:3217), rear_spar(1:3217), 'r-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(-x_te(1:3101), front_spar(1:3101), 'r-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(-x_te(1:3217), rear_spar(1:3217), 'r-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot([-9.75 , 9.75] , [rear_spar(976) , rear_spar(976)] , 'r-', 'LineWidth', 1, 'DisplayName' , 'Secondary Spar');

% Plotting heavy ribs
plot(rib_1(1,:) , rib_1(2,:) , 'g-' , 'LineWidth', 1.5 , 'DisplayName', 'Heavy Ribs')
plot(rib_2(1,:) , rib_2(2,:) , 'g-' , 'LineWidth', 1.5 , 'HandleVisibility', 'off')
plot(rib_3(1,:) , rib_3(2,:) , 'g-' , 'LineWidth', 1.5 , 'HandleVisibility', 'off')
plot(rib_4(1,:) , rib_4(2,:) , 'g-' , 'LineWidth', 1.5 , 'HandleVisibility', 'off')
plot(rib_5(1,:) , rib_5(2,:) , 'g-' , 'LineWidth', 1.5 , 'HandleVisibility', 'off')
plot(rib_6(1,:) , rib_6(2,:) , 'g-' , 'LineWidth', 1.5 , 'HandleVisibility', 'off')
plot(-rib_1(1,:) , rib_1(2,:) , 'g-' , 'LineWidth', 1.5 , 'HandleVisibility', 'off')
plot(-rib_2(1,:) , rib_2(2,:) , 'g-' , 'LineWidth', 1.5 , 'HandleVisibility', 'off')
plot(-rib_3(1,:) , rib_3(2,:) , 'g-' , 'LineWidth', 1.5 , 'HandleVisibility', 'off')
plot(-rib_4(1,:) , rib_4(2,:) , 'g-' , 'LineWidth', 1.5 , 'HandleVisibility', 'off')
plot(-rib_5(1,:) , rib_5(2,:) , 'g-' , 'LineWidth', 1.5 , 'HandleVisibility', 'off')
plot(-rib_6(1,:) , rib_6(2,:) , 'g-' , 'LineWidth', 1.5 , 'HandleVisibility', 'off')

% Filling the flaps (right side) with transparency
fill(inner_flap(1,:), inner_flap(2,:), 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'DisplayName', 'Flaps'); % Transparent cyan flap
fill(outer_flap(1,:), outer_flap(2,:), 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'HandleVisibility', 'off'); % Transparent cyan flap

% Filling the flaps (left side) with transparency
fill(-inner_flap(1,:), inner_flap(2,:), 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'HandleVisibility', 'off'); % Transparent cyan flap
fill(-outer_flap(1,:), outer_flap(2,:), 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'HandleVisibility', 'off'); % Transparent cyan flap

%filling the ailerons
fill(AILERON(1,:), AILERON(2,:), [0.9290 0.6940 0.1250] , 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'DisplayName','Ailerons'); % Transparent cyan flap
fill(-AILERON(1,:), AILERON(2,:), [0.9290 0.6940 0.1250] , 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'HandleVisibility', 'off'); % Transparent cyan flap


%adding slats
fill(inner_slat(1,:), inner_slat(2,:), 'y', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'DisplayName', 'Slats'); % Transparent cyan flap
fill(centre_slat(1,:), centre_slat(2,:), 'y', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'HandleVisibility', 'off'); % Transparent cyan flap
fill(outer_slat(1,:), outer_slat(2,:), 'y', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'HandleVisibility', 'off'); % Transparent cyan flap

% Filling the flaps (left side) with transparency
fill(-inner_slat(1,:), inner_slat(2,:), 'y', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'HandleVisibility', 'off'); % Transparent cyan flap
fill(-centre_slat(1,:), centre_slat(2,:), 'y', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'HandleVisibility', 'off'); % Transparent cyan flap
fill(-outer_slat(1,:), outer_slat(2,:), 'y', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'HandleVisibility', 'off'); % Transparent cyan flap

%plotting tanks
fill(inner_tank(1,:), inner_tank(2,:), 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'DisplayName', 'Fuel Tanks');
fill(centre_tank(1,:), centre_tank(2,:), 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'HandleVisibility', 'off');
fill(outer_tank(1,:), outer_tank(2,:), 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'HandleVisibility', 'off');

fill(-inner_tank(1,:), inner_tank(2,:), 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'HandleVisibility', 'off');
fill(-centre_tank(1,:), centre_tank(2,:), 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'HandleVisibility', 'off');
fill(-outer_tank(1,:), outer_tank(2,:), 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'HandleVisibility', 'off');

%central tank
%fill(fuselage_tank(1,:), fuselage_tank(2,:), [0.2 0.2 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'DisplayName' , 'Fuselage Tank');

%plotting undercarriage
plot(x_values , y_values , 'r--', 'HandleVisibility', 'off')
fill(x_values, y_values , 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'DisplayName','Undercarriage Retraction')
plot(-x_values , y_values , 'r--', 'HandleVisibility', 'off')
fill(-x_values, y_values , 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'HandleVisibility', 'off')

%plotting fuselage
plot(fuselage(1,:) , fuselage(2,:) , 'Color' , [0.2 0.2 0.2], 'LineStyle', '--' , 'LineWidth' , 0.5 , 'DisplayName' , 'Fuselage (Max Width)')
plot([2.55 , 2.55] , fuselage(2,:) , 'Color' , [0.2 0.2 0.2], 'LineStyle', ':' , 'LineWidth' , 0.5 , 'DisplayName' , 'Fuselage (Wing Fuselage Interface)')
plot(-fuselage(1,:) , fuselage(2,:) , 'Color' , [0.2 0.2 0.2] , 'LineStyle',  '--' , 'LineWidth' , 0.5 ,  'HandleVisibility', 'off')
plot([-2.55 , -2.55] , fuselage(2,:) , 'Color' , [0.2 0.2 0.2] , 'LineStyle', ':' , 'LineWidth' , 0.5 ,  'HandleVisibility', 'off')


% Remove ticks but keep the grid
grid on;
set(gca, 'XTickLabel', [], 'YTickLabel', []);

% Adjust axis limits and aspect ratio
xlim([-35 50]); % Adjust to allow legend space on the right
ylim([-20 10]); % Keep consistent vertical limits
pbaspect([3 1 1]); % Aspect ratio to ensure appropriate width for the legend

% Legend
legend('FontSize', 15);

% Box around the plot
box on;
hold off;




%{
% draw the wing
figure
hold on
plot([0, -wing_span / 2], [le_root, le_tip], 'b-', 'LineWidth', 1)
plot([0, -kink_span], [te_root, te_root], 'b-', 'LineWidth', 1)
plot([-kink_span, -wing_span / 2], [te_root, te_tip], 'b-', 'LineWidth', 1)
plot([-wing_span / 2, -wing_span / 2], [le_tip, te_tip], 'b-', 'LineWidth', 1)
plot([0, wing_span / 2], [le_root, le_tip], 'b-', 'LineWidth', 1)
plot([0, kink_span], [te_root, te_root], 'b-', 'LineWidth', 1)
plot([kink_span, wing_span / 2], [te_root, te_tip], 'b-', 'LineWidth', 1)
plot([wing_span / 2, wing_span / 2], [le_tip, te_tip], 'b-', 'LineWidth', 1)

grid on
axis equal
%}
%%

clear
%Tailplane planform
%horizontal

% parameters
area = 29.08 * 2; % m^2
AR = 1.6;
wing_span = (area * AR * 2) ^ 0.5; % m
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

grad_le = (le_tip - le_root)/(wing_span/2);
grad_te = (te_tip - te_root)/(wing_span/2);

x_vals_le = 1.56:0.01:wing_span/2;

y_le = grad_le * x_vals_le + le_root;

x_vals_te = 0.61:0.01:wing_span/2;

y_te = grad_te * x_vals_te + te_root;

abox = [-0.61 , 0.61 , 1.56 , -1.56 , -0.61 ; y_te(1) , y_te(1) , y_le(1) , y_le(1) , y_te(1)];

front_spar = [0 , wing_span/2 ; le_root - (0.25 * root_chord) , le_tip - (0.25 * tip_chord)];

chord_grad = (tip_chord-root_chord)/(wing_span/2);

x_fs = 1.32:0.01:wing_span/2;

chord_dist = x_fs * chord_grad + root_chord;

front_spar = y_te(72:end) + 0.75 * chord_dist;

front_spar = [front_spar(1) , front_spar];

x_fs = [0 , x_fs];

rear_spar = [0 , wing_span/2 ; le_root - (0.65 * root_chord) , le_tip - (0.65 * tip_chord)];

x_rs = 0.94:0.01:wing_span/2;

chord_dist = x_rs * chord_grad + root_chord;

rear_spar = y_te(34:end) + 0.35 * chord_dist;

rear_spar = [rear_spar(1) , rear_spar];

x_rs = [0 , x_rs];

% draw the wing
figure
clf
hold on

%{
plot([0, wing_span / 2], [le_root, le_tip], 'k-', 'LineWidth', 1)
plot([0, wing_span / 2], [te_root, te_tip], 'k-', 'LineWidth', 1)
plot([wing_span / 2, wing_span / 2], [le_tip, te_tip], 'k-', 'LineWidth', 1)
plot(-[0, wing_span / 2], [0, qcy], 'r--', 'LineWidth', 1)
plot(-[0, wing_span / 2], [le_root, le_tip], 'k-', 'LineWidth', 1)
plot(-[0, wing_span / 2], [te_root, te_tip], 'k-', 'LineWidth', 1)
plot(-[wing_span / 2, wing_span / 2], [le_tip, te_tip], 'k-', 'LineWidth', 1)
%}

plot(x_vals_le , y_le , 'k-', 'LineWidth', 1 , 'DisplayName', 'Horizontal Stabiliser Outline')
plot(x_vals_te , y_te , 'k-', 'LineWidth', 1 , 'HandleVisibility', 'off')
plot(-x_vals_le , y_le , 'k-', 'LineWidth', 1 , 'HandleVisibility', 'off')
plot(-x_vals_te , y_te , 'k-', 'LineWidth', 1 , 'HandleVisibility', 'off')
plot([wing_span / 2, wing_span / 2], [le_tip, te_tip], 'k-', 'LineWidth', 1 , 'HandleVisibility', 'off')
plot(-[wing_span / 2, wing_span / 2], [le_tip, te_tip], 'k-', 'LineWidth', 1 , 'HandleVisibility', 'off')

plot(abox(1,:),abox(2,:), 'Color' , [0.9290 0.6940 0.1250] , 'LineStyle' , '--' , 'DisplayName', 'Fuselage')

plot(x_fs,front_spar , 'r' , 'LineWidth',2 , 'DisplayName', 'Spars')
plot(-x_fs,front_spar , 'r' , 'LineWidth',2 , 'HandleVisibility', 'off')
plot(x_rs,rear_spar , 'r' , 'LineWidth',2 , 'HandleVisibility', 'off')
plot(-x_rs,rear_spar , 'r' , 'LineWidth', 2 , 'HandleVisibility', 'off')

% Remove ticks but keep the grid
grid on;
set(gca, 'XTickLabel', [], 'YTickLabel', []);

axis equal
xlim([-7.5 7.5])
ylim([-10 5])

legend(FontSize=15);
% Box around the plot
box on;
hold off;

