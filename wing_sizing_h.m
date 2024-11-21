clear
clf
clc

% parameters
area = 68; % m^2
AR = 5.8;
wing_span = (area * AR / 2) ^ 0.5; % m
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

