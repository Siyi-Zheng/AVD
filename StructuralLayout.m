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

%inner flaps


%spar readjustment (not needed now)
%where kink spar and central box spar meet
grad_rspar = (y_le(end) - y_le(976)) / (x_le(end) - 9.75);
c = rear_spar(976) - grad_rspar * 9.75;
intersect = (front_spar(1) - c) / grad_rspar;

%spar le angle
spar_fangle = atand( ( front_spar(256) - front_spar(end) ) / ( x_le(end) - 2.55) );
spar_flength = (x_le(end) - 2.55) / cosd(spar_fangle);

%rear spar calculations
dist = (17.9543 + (y_le(1) - rear_spar(1))) * 1000;
spar_rangle1 = atand( ( rear_spar(256) - rear_spar(976) ) / ( x_le(976) - 2.55) );
spar_rlength1 = (x_le(976) - 2.55) / cosd(spar_rangle1);
spar_rangle2 = atand( ( rear_spar(976) - rear_spar(end) ) / ( x_le(end) - 9.75) );
spar_rlength2 = (x_le(end) - 9.75) / cosd(spar_rangle2);

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
%plotting heavy ribs
%rib one at 2.55
rib_1 = [2.55 , 2.55 ; y_le(256) , y_te(256)];


%rib 2 at 30% span
find = 0.3 * wing_span/2;
%index = find(x_le == find);
rib_2 = [find , find ; y_le(976) , y_te(976)];



%rib 3 at 65% span
find = round(0.65 * wing_span/2 , 2);
index = 2114;
rib_3 = [find , find ; y_le(index) , y_te(index)];

%rib 4 at end of spars
rib_4 = [x_le(end) , x_le(end) ; y_le(end) , y_te(2973)];



%%
%Flaps
%front line of flaps is 25% chord, and flaps range from 10% span to 65%
span_65 = 0.65 * wing_span/2;
span_10 = 0.1 * wing_span/2;

index1 = round(span_10 * 100) + 2;
index2 = round(span_65 * 100) - 1;

flap_le = y_le(index1:index2) - 0.75 .* chord(index1:index2);
flap_te = y_te(index1:index2);

FLAP = [ x_le(index1:index2) , fliplr(x_le(index1:index2)) , x_le(index1) ; flap_le , fliplr(flap_te) , flap_le(1) ];

%split flap to inner and outer
inner_flap = [FLAP(:,1:635) , FLAP(:, 2938:end)];
outer_flap = FLAP(:,665:2908);

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



%% PLOTTING

figure
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
plot(x_le, front_spar, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Spar'); % Included in legend
plot(x_le, rear_spar, 'r-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(-x_le, front_spar, 'r-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(-x_le, rear_spar, 'r-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot([-9.75 , 9.75] , [rear_spar(976) , rear_spar(976)] , 'r-', 'LineWidth', 1, 'HandleVisibility', 'off');

% Plotting heavy ribs
plot(rib_1(1,:) , rib_1(2,:) , 'g-' , 'LineWidth', 1.5 , 'DisplayName', 'Heavy Ribs')
plot(rib_2(1,:) , rib_2(2,:) , 'g-' , 'LineWidth', 1.5 , 'HandleVisibility', 'off')
plot(rib_3(1,:) , rib_3(2,:) , 'g-' , 'LineWidth', 1.5 , 'HandleVisibility', 'off')
plot(rib_4(1,:) , rib_4(2,:) , 'g-' , 'LineWidth', 1.5 , 'HandleVisibility', 'off')
plot(-rib_1(1,:) , rib_1(2,:) , 'g-' , 'LineWidth', 1.5 , 'HandleVisibility', 'off')
plot(-rib_2(1,:) , rib_2(2,:) , 'g-' , 'LineWidth', 1.5 , 'HandleVisibility', 'off')
plot(-rib_3(1,:) , rib_3(2,:) , 'g-' , 'LineWidth', 1.5 , 'HandleVisibility', 'off')
plot(-rib_4(1,:) , rib_4(2,:) , 'g-' , 'LineWidth', 1.5 , 'HandleVisibility', 'off')

% Filling the flaps (right side) with transparency
fill(inner_flap(1,:), inner_flap(2,:), 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'DisplayName', 'Flaps'); % Transparent cyan flap
fill(outer_flap(1,:), outer_flap(2,:), 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'HandleVisibility', 'off'); % Transparent cyan flap

% Filling the flaps (left side) with transparency
fill(-inner_flap(1,:), inner_flap(2,:), 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'HandleVisibility', 'off'); % Transparent cyan flap
fill(-outer_flap(1,:), outer_flap(2,:), 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'HandleVisibility', 'off'); % Transparent cyan flap


%adding slats
fill(inner_slat(1,:), inner_slat(2,:), 'y', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'DisplayName', 'Slats'); % Transparent cyan flap
fill(centre_slat(1,:), centre_slat(2,:), 'y', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'HandleVisibility', 'off'); % Transparent cyan flap
fill(outer_slat(1,:), outer_slat(2,:), 'y', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'HandleVisibility', 'off'); % Transparent cyan flap

% Filling the flaps (left side) with transparency
fill(-inner_slat(1,:), inner_slat(2,:), 'y', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'HandleVisibility', 'off'); % Transparent cyan flap
fill(-centre_slat(1,:), centre_slat(2,:), 'y', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'HandleVisibility', 'off'); % Transparent cyan flap
fill(-outer_slat(1,:), outer_slat(2,:), 'y', 'EdgeColor', 'none', 'FaceAlpha', 0.5 , 'HandleVisibility', 'off'); % Transparent cyan flap

grid on
axis equal
legend;
hold off




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



