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

grad_le = (le_tip-le_root)/((wing_span)/2);
grad_te = (te_tip-te_root)/(((wing_span)/2)-kink_span);

x_r = 0:0.01:(wing_span/2);
x_l = -x_r;

y_le = x_r * grad_le + le_root;
y_te = (x_r(1:find(x_r == 9.75)) * 0) + te_root;
c = te_root - grad_te * kink_span;
y_te = ([y_te , (x_r(find(x_r == 9.76):end)) * grad_te + c]);

y_le_box = 

%adding central box
x_box = -2.55:0.01:2.55;
y_lebox = (x_box * 0) + 2.55 * grad_le + le_root;
y_tebox = (x_box * 0) + te_root;




%adding spars
%le_spar @ 12% chord
chord = y_le - y_te;
front_spar = y_le - 0.12 * chord;
rear_spar = y_le - 0.7 * chord;


figure
hold on

%plotting wing
plot(x_r , y_le , 'k-' , LineWidth=1.5)
plot(x_r , y_te , 'k-' , LineWidth=1.5)
plot(x_l , y_le , 'k-' , LineWidth=1.5)
plot(x_l , y_te , 'k-' ,  LineWidth=1.5)
plot([-wing_span / 2, -wing_span / 2], [le_tip, te_tip], 'k-', 'LineWidth', 1.5)
plot([wing_span / 2, wing_span / 2], [le_tip, te_tip], 'k-', 'LineWidth', 1.5)

%plotting wing box
plot(x_box , y_lebox , 'Color' , [0.9290 0.6940 0.1250], 'LineStyle', '--', 'LineWidth', 1.5)
plot(x_box , y_tebox , 'Color' , [0.9290 0.6940 0.1250], 'LineStyle', '--', 'LineWidth', 1.5)
plot([-2.55 , -2.55] , [y_lebox(1) , y_tebox(1)] , 'Color' , [0.9290 0.6940 0.1250] , 'LineStyle', '--', 'LineWidth', 1.5)
plot([2.55 , 2.55] , [y_lebox(end) , y_tebox(end)] , 'Color' , [0.9290 0.6940 0.1250] , 'LineStyle', '--', 'LineWidth', 1.5)

%plotting spars
plot(x_r , front_spar , 'r-' , LineWidth=1.5)
plot(x_r , rear_spar , 'r-' , LineWidth=1.5)
plot(-x_r , front_spar , 'r-' , LineWidth=1.5)
plot(-x_r , rear_spar , 'r-' , LineWidth=1.5)

grid on
axis equal
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



