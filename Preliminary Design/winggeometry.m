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

%adding central box
x_box = -2.55:0.01:2.55;
y_leboxa = (x_box * 0) + 2.55 * grad_le + le_root;
y_teboxa = (x_box * 0) + te_root;


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
plot(x_box, y_leboxa, 'Color', 'k' , 'LineStyle', '-', 'LineWidth', 1.5, 'DisplayName', 'WingBox'); % Included in legend
plot(x_box, y_teboxa, 'Color', 'k' , 'LineStyle', '-', 'LineWidth', 1.5, 'HandleVisibility', 'off');

axis equal
hold off