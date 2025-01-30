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
plot([0, wing_span / 2], [le_root, le_tip], 'k-', 'LineWidth', 1)
plot([0, kink_span], [te_root, te_root], 'k-', 'LineWidth', 1)
plot([kink_span, wing_span / 2], [te_root, te_tip], 'k-', 'LineWidth', 1)
plot([0, 0], [le_root, te_root], 'k-', 'LineWidth', 1)
plot([wing_span / 2, wing_span / 2], [le_tip, te_tip], 'k-', 'LineWidth', 1)
axis equal

% get chord and spars
x = linspace(0, wing_span/2, 1000);
le = le_root + (le_tip - le_root) .* x ./ (wing_span / 2);
te1 = te_root + (te_root - te_tip) .* (kink_span - x) ./ ((wing_span / 2) - kink_span);
te = min(te_root, te1);
chord = le - te;
% front spar/ 
front_spar = le - 0.12 * chord;
rear_spar = te + 0.3 * chord;
front_spar = linspace(front_spar(1), front_spar(end), 1000);
plot(x, front_spar, "r-")
plot(x, rear_spar, "r-")

% stringer layout
tip_stringers = 16; % number of stringers at the top
stringer_spacing = (front_spar(end) - rear_spar(end)) / (tip_stringers + 1);
for i = 1:tip_stringers * 6
    stringer_pos = max(rear_spar, front_spar - i * stringer_spacing);
    end_pos = 1000;
    if any(stringer_pos == rear_spar)
        end_pos = find(stringer_pos == rear_spar, 1);
    end
    plot(x(1:end_pos), stringer_pos(1:end_pos), "b-")
end

% get the number of stringers at various locations
num_stringers = floor((front_spar - rear_spar) / stringer_spacing);
num_stringers_root = num_stringers(1);
num_stringers_kink = num_stringers(ceil(kink_span ./ wing_span * 2000));
