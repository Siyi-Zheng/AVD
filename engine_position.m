span = 65/2; 
engine_mass = 6147;
engine_thrust = 341200; % Thrust per engine in Newtons
mean_chord = 7.41; % Mean chord length in meters 
fan_diameter = 266; 
wing_mass_per_unit_span = 500; % Wing mass per meter of span in kg, average of a wide body
total_lift = (engine_mass * 9.81 * 4); 

% Discretize the wing span
num_points = 200;
x_span = linspace(0, span/2, num_points); % Only half-span needed due to symmetry

% Define possible engine position ranges
engine_positions = linspace(0.1 * (span / 2), 0.9 * (span / 2), 50); % Range from 10% to 90% of half-span


q_aero = 240980 * sqrt(1 - (x_span / (span / 2)).^2); % Aerodynamic load inspired by document
q_wing = 318825 * (x_span / (span / 2)); % Wing weight distribution
q_fuel = 96392 * (x_span / (span / 2)); % Fuel weight distribution

% Total distributed load
q_total = q_aero + q_wing + q_fuel;

% Calculate Bending Moment due to Total Distributed Load
bending_moment = cumtrapz(x_span, q_total .* (x_span(end) - x_span));

min_bending_moment = inf; 
optimal_positions = [0, 0]; 


optimal_inertia = inf; 
optimal_yaw_moment = inf; 

for pos1 = engine_positions
    for pos2 = engine_positions
        % Ensure pos2 is at least 10 meters after pos1
        if pos2 > pos1 + 10
            % Calculate moment contributions from engines (weight and thrust)
            weight_moment_engine1 = engine_mass * 9.81 * (span / 2 - pos1); 
            weight_moment_engine2 = engine_mass * 9.81 * (span / 2 - pos2);
            thrust_moment_engine1 = engine_thrust * (span / 2 - pos1);
            thrust_moment_engine2 = engine_thrust * (span / 2 - pos2);
            
            % Calculate total bending moment at the root
            total_moment = bending_moment + weight_moment_engine1 + weight_moment_engine2 + thrust_moment_engine1 + thrust_moment_engine2;
            root_bending_moment = total_moment(1); 
            
            % Calculate moment of inertia for current engine positions
            inertia_engine1 = engine_mass * (span / 2 - pos1)^2;
            inertia_engine2 = engine_mass * (span / 2 - pos2)^2;
            total_inertia = inertia_engine1 + inertia_engine2;
            
            % Calculate yaw moment due to engine failure
            % Assume failure of engine at pos2 for worst case
            yaw_moment_failure = engine_thrust * (span / 2 - pos2);
            
            
            combined_score = root_bending_moment + 0.1 * total_inertia + 0.05 * yaw_moment_failure;
            
            % Update optimal positions based on combined score
            if combined_score < min_bending_moment
                min_bending_moment = combined_score;
                optimal_positions = [pos1, pos2];
                optimal_inertia = total_inertia;
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
fprintf('Minimum Bending Moment: %.2f Nm\n', min_bending_moment);
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
plot(optimal_positions(1), min_bending_moment, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
plot(optimal_positions(2), min_bending_moment, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
xlabel('Engine Position (m)');
ylabel('Total Bending Moment at Root (Nm)');
title('Total Bending Moment vs. Engine Position');
legend('Bending Moment', 'Optimal Engine Position 1', 'Optimal Engine Position 2');
