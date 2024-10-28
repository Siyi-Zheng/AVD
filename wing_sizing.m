area = 481.77; % m^2
wing_span = 65; % m
taper_ratio = 0.25;
quarter_chord_sweep = 30; % degrees
trailing_edge_kink = 0.2; % fraction of wing span

% simple design
total_chord = 1 + 1 / taper_ratio;
tip_chord = 2 * area / (wing_span * total_chord);
root_chord = tip_chord / taper_ratio;
qcy = -wing_span * tand(quarter_chord_sweep);
le_tip = qcy + tip_chord / 4;
te_tip = qcy - 3 * tip_chord / 4;
le_root = root_chord / 4;
te_root = -3 * root_chord / 4;
% draw the wing
figure
clf
hold on
plot([0, wing_span], [0, qcy], 'r--', 'LineWidth', 1)
plot([0, wing_span], [le_root, le_tip], 'k-', 'LineWidth', 1)
plot([0, wing_span], [te_root, te_tip], 'k-', 'LineWidth', 1)
plot([0, 0], [le_root, te_root], 'k-', 'LineWidth', 1)
plot([wing_span, wing_span], [le_tip, te_tip], 'k-', 'LineWidth', 1)

% wing with trailing edge kink
kink_span = wing_span * trailing_edge_kink;
extra_area = 0.5 * kink_span ^ 2 * tand(quarter_chord_sweep);
new_area = area - extra_area;
% recalculate the wing
tip_chord2 = 2 * new_area / (wing_span * total_chord);
root_chord2 = tip_chord2 / taper_ratio + kink_span * tand(quarter_chord_sweep);
qcy = -wing_span * tand(quarter_chord_sweep);
le_tip = qcy + tip_chord2 / 4;
te_tip = qcy - 3 * tip_chord2 / 4;
le_root = root_chord2 / 4;
te_root = -3 * root_chord2 / 4;
% draw the wing
hold on
plot([0, -wing_span], [0, qcy], 'r--', 'LineWidth', 1)
plot([0, -wing_span], [le_root, le_tip], 'b-', 'LineWidth', 1)
plot([0, -kink_span], [te_root, te_root], 'b-', 'LineWidth', 1)
plot([-kink_span, -wing_span], [te_root, te_tip], 'b-', 'LineWidth', 1)
plot([0, 0], [le_root, te_root], 'b-', 'LineWidth', 1)
plot([-wing_span, -wing_span], [le_tip, te_tip], 'b-', 'LineWidth', 1)
axis equal
