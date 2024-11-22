clear
clf
clc

% parameters
area = 482; % m^2
wing_span = 65; % m
taper_ratio = 0.25;
quarter_chord_sweep = 26.6; % degrees
trailing_edge_kink = 0.3; % fraction of wing span
area_ratio = 0.0939; % cross-sectional area as a proportion of c^2

% simple design
total_chord = 1 + 1 / taper_ratio;
tip_chord = 2 * area / (wing_span * total_chord);
root_chord = tip_chord / taper_ratio;
qcy = -wing_span / 2 * tand(quarter_chord_sweep);
le_tip = qcy + tip_chord / 4;
te_tip = qcy - 3 * tip_chord / 4;
le_root = root_chord / 4;
te_root = -3 * root_chord / 4;

% draw the wing
figure
clf
hold on
plot([0, wing_span / 2], [0, qcy], 'r--', 'LineWidth', 1)
plot([0, wing_span / 2], [le_root, le_tip], 'k-', 'LineWidth', 1)
plot([0, wing_span / 2], [te_root, te_tip], 'k-', 'LineWidth', 1)
plot([0, 0], [le_root, te_root], 'k-', 'LineWidth', 1)
plot([wing_span / 2, wing_span / 2], [le_tip, te_tip], 'k-', 'LineWidth', 1)
axis equal

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

% draw the wing
hold on
plot([0, -wing_span / 2], [0, qcy], 'r--', 'LineWidth', 1)
plot([0, -wing_span / 2], [le_root, le_tip], 'b-', 'LineWidth', 1)
plot([0, -kink_span], [te_root, te_root], 'b-', 'LineWidth', 1)
plot([-kink_span, -wing_span / 2], [te_root, te_tip], 'b-', 'LineWidth', 1)
plot([0, 0], [le_root, te_root], 'b-', 'LineWidth', 1)
plot([-wing_span / 2, -wing_span / 2], [le_tip, te_tip], 'b-', 'LineWidth', 1)
axis equal

figure
x = linspace(3.17, wing_span/2, 1000);
% get leading edge coordinates
le = le_root + (le_tip - le_root) .* x ./ (wing_span / 2);
te1 = te_root + (te_root - te_tip) .* (kink_span - x) ./ ((wing_span / 2) - kink_span);
te = min(te_root, te1);
chord = le - te;
csa = chord .^ 2 * area_ratio;
plot(x, csa)
wing_vol = trapz(x, csa);
