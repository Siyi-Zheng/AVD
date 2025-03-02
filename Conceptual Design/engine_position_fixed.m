clear
clc

span = 32.5; % actually semi-span due to symmetry (m)
engine_mass = 9000; % (kg) - old value 6147 but need to take nacelle mass into account
engine_thrust = 299800; % Thrust per engine (N)
mean_chord = 7.41; % Mean chord length in meters 
fan_diameter = 2.66; % (m), including nacelle
lift = 390000 * 9.81; % (N) also the weight of the plane
wing_weight = 43000 * 9.81; % (N) for a single wing - just guesses for now
single_fuel_weight = 90000 * 9.81;  % (N) fuel in a single wing - just guesses for now

% Discretize the wing span
num_points = 200;
x_span = linspace(0, span, num_points); % Only semi-span needed due to symmetry

% Define possible engine position ranges
engine_positions = linspace(0.25 * (span), 0.75 * (span), 200); % Range of sensible values

q_aero = lift / 2 * sqrt(1 - (x_span / (span)).^2) / (pi * span / 4); % Aerodynamic load inspired by document
q_wing = wing_weight * ((span - x_span) / (span)) / (span / 2); % Wing weight distribution
q_fuel = single_fuel_weight * ((span - x_span) / (span)) / (span / 2); % Fuel weight distribution
wing_mass_per_unit_span = (q_wing + q_fuel) / 9.81;
inertia_wing = trapz(x_span, wing_mass_per_unit_span .* x_span.^2); % moment of inertia of wing

% Total distributed load
q_total = q_aero - q_wing - q_fuel;

% Calculate Bending Moment due to Total Distributed Load
bending_moment = x_span; % preallocate to the same size
for point = 1:length(x_span)
    % integrate all points to the right of the point
    moment = sum(x_span(point:end) .* q_total(point:end)) * (span /( length(x_span) - 1));
    bending_moment(point) = moment;
end

optimal_score = inf;
optimal_positions = [0, 0]; 

optimal_moment = inf;
optimal_inertia = inf; 
optimal_yaw_moment = inf; 

for pos1 = engine_positions
    for pos2 = engine_positions
        % Ensure pos2 is at least 10 meters after pos1
        if pos2 > pos1 + 10
            % Calculate moment contributions from engines (weight and thrust)
            weight_moment_engine1 = -engine_mass * 9.81 * (span - pos1); 
            weight_moment_engine2 = -engine_mass * 9.81 * (span - pos2);
            
            % Calculate total bending moment at the root
            total_moment = bending_moment + weight_moment_engine1 + weight_moment_engine2;
            root_bending_moment = total_moment(1); 
            
            % Calculate moment of inertia for current engine positions
            inertia_engine1 = engine_mass * (span - pos1)^2;
            inertia_engine2 = engine_mass * (span - pos2)^2;
            total_inertia = inertia_engine1 + inertia_engine2 + inertia_wing;
            
            % Calculate yaw moment due to engine failure
            % Assume failure of engine at pos2 for worst case
            yaw_moment_failure = engine_thrust * (span - pos2);
            
            combined_score = root_bending_moment / optimal_moment + ...
                0.2 * total_inertia / optimal_inertia + ...
                0.3 * yaw_moment_failure / optimal_yaw_moment;
            if combined_score == 0
                combined_score = 3; % to fix the first iteration
            end
            
            % Update optimal positions based on combined score
            if combined_score < optimal_score
                optimal_score = combined_score;
                optimal_positions = [pos1, pos2];
            end
            if root_bending_moment < optimal_moment
                optimal_moment = root_bending_moment;
            end
            if total_inertia < optimal_inertia
                optimal_inertia = total_inertia;
            end
            if yaw_moment_failure < optimal_yaw_moment
                optimal_yaw_moment = yaw_moment_failure;
            end
        end
    end
end

% Account for forward and downward placement of engines relative to leading edge
forward_position = 2 * fan_diameter; % 2 fan diameters forward from the leading edge
downward_position = -1 * fan_diameter; % 1 fan diameter down from the leading edge

% Display Optimal Positions and Calculated Values
fprintf('Optimal Engine Positions along the Wing Span: %.2f meters and %.2f meters\n', optimal_positions(1), optimal_positions(2));
fprintf('Forward Position (relative to leading edge): %.2f meters\n', forward_position);
fprintf('Downward Position (relative to leading edge): %.2f meters\n', downward_position);
fprintf('Minimum Bending Moment: %.2f Nm\n', optimal_moment);
fprintf('Optimal Moment of Inertia: %.2f kg*m^2\n', optimal_inertia);
fprintf('Yaw Moment due to Engine Failure: %.2f Nm\n', optimal_yaw_moment);

% Plot Distributed Loads
figure;
subplot(3,1,1);
plot(x_span, q_aero, 'b', 'LineWidth', 1.5);
xlabel('Spanwise Position (m)');
ylabel('Aerodynamic Load (N/m)');
title('Distributed Aerodynamic Load');

subplot(3,1,2);
plot(x_span, q_wing + q_fuel, 'r', 'LineWidth', 1.5);
xlabel('Spanwise Position (m)');
ylabel('Structural + Fuel Load (N/m)');
title('Distributed Load from Wing and Fuel');

subplot(3,1,3);
plot(x_span, q_total, 'k', 'LineWidth', 1.5);
xlabel('Spanwise Position (m)');
ylabel('Total Load (N/m)');
title('Total Distributed Load Along Wingspan');

% Plot Results with Optimal Positions
figure;
subplot(2, 1, 1);
plot(x_span, q_total, 'b', 'LineWidth', 1.5);
xlabel('Spanwise Position (m)');
ylabel('Total Distributed Load (N/m)');
title('Total Distributed Load along the Wing Span');

subplot(2, 1, 2);
plot(engine_positions, bending_moment, 'r', 'LineWidth', 1.5);
hold on;
plot(optimal_positions(1), optimal_moment, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
plot(optimal_positions(2), optimal_moment, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
xlabel('Engine Position (m)');
ylabel('Total Bending Moment at Root (Nm)');
title('Total Bending Moment vs. Engine Position');
legend('Bending Moment', 'Optimal Engine Position 1', 'Optimal Engine Position 2');
