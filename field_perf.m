% Constants and Aircraft Parameters
rho = 1.225  %Air density at sea level (kg/m^3)
S_ref = 482;
W_S = W/S_ref;  %Wing loading (N/m^2)
mu = 0.03  %Rolling friction coefficient
CL = 2.15  %Lift coefficient
CD_0 = 0.1026  %Zero-lift drag coefficient
AR = 8.77  %Aspect ratio
e = 0.95  %Oswald efficiency factor
T = 1123e3  %Thrust (N)
V_S = 70.621  %Stall speed (m/s)
W = 390000*9.81  %Aircraft weight (N)
W_landing = 0.85 * W  %Landing weight (N), 85% of MTOW
V1 = 0  %Inititial velocity (m/s)
V2 = 1.1 * V_S  %Final velocity for ground roll (m/s)
L_D = 0.2  %Lift-to-drag ratio

% Takeoff Ground Roll Distance Calculation
K_A = (rho / (2 * W_S)) * (mu * CL - CD_0 - (CL^2) / (pi * AR * e));
K_T = T / W - mu;

S_g = (1 / (2 * 9.81 * K_A)) * log((K_T + K_A * V2^2) / (K_T + K_A * V1^2));

% Rotation Distance
alpha_lof = 12;       % Angle of attack at lift-off (degrees)
alpha_G = 0;            % Ground angle of attack (degrees)
dtheta_dt = 3;          % Pitch rate (degrees/s)
V_lof = V2;             % Lift-off velocity (m/s)

S_R = ((alpha_lof - alpha_G) / dtheta_dt) * V_lof;

% Transition and Climb Distance

V_TR = 1.15 * V_S;       % Transition velocity
n = 3.75;                 % Load factor
R = (V_TR^2) / ((n - 1) * 9.81); % Transition radius (m)
gamma_CL = asin(T / W - 1 / L_D); % Climb angle (radians)

h_TR = R * (1 - cos(gamma_CL)); % Transition height (m)
h_OBS = 10.67;                      % Obstacle height (m)

if h_OBS >= h_TR
    S_TR = R * sin(gamma_CL);
    S_CL = (h_OBS - h_TR) / tan(gamma_CL);
else
    S_TR = sqrt(R^2 - (R - h_OBS)^2);
    S_CL = 0;
end

% Total Takeoff Distance
S_TO = 1.15 * (S_g + S_R + S_TR + S_CL);

% Balanced Field Length (BFL)
N_e = 4;                   % Number of engines
T_W_O = T / W;             % Thrust-to-weight ratio
D_OEI = 2e4;               % Drag with one engine inoperative (N)
CL_max_clean = 1.51;       % Max lift coefficient in clean configuration
BPR = 9.6;                 % Bypass ratio
T_takeoff = T;
CL_climb = 1;              % Lift coefficient during climb

sigma = 1;                 % Air density ratio
T_dash = 0.75 * T_takeoff * ((5 + BPR) / (4 + BPR));
U = 0.01 * CL_max_clean + 0.02;
gamma_min = 0.018 + 0.003 * N_e;
gamma_CL = asin(((N_e - 1) / N_e) * T_W_O - D_OEI / W);
G = gamma_CL - gamma_min;

BFL = 0.863 / (1 + 2.3 * G) * (W_S / (rho * 9.81 * CL_climb + h_OBS) * (1 / ((T_dash / W) - U) + 2.7)) + 655 / sqrt(sigma);

%Landing Parameters with Updated Obstacle Height (50 ft = 15.24 m)
CL_max_landing = 2.85 
h_OBS_landing = 15.24  %Obstacle height for landing (m)
CD_0 = 0.1424  %Zero-lift drag coefficient during landing
V_S_landing = sqrt(W_landing/(rho*S_ref*CL_max_landing))  %Adjusted stall speed for landing weight
V_a = 1.3 * V_S_landing  %Approach speed (m/s)
gamma_a = deg2rad(2)  %Approach angle (radians)
V_F = 1.23 * V_S_landing  %Final approach speed (m/s)
R = (V_F^2) / ((n - 1) * 9.81)  %Radius for flare phase (m)
h_f = R * (1 - cos(gamma_a))  %Flare height (m)

S_a = (h_OBS_landing - h_f) / tan(gamma_a)
S_F = R * sin(gamma_a)

%Free Roll and Braking Distance
V_TD = 1.15 * V_S %V_S_landing  # Touchdown speed (m/s)
t_FR = 2  %Free roll time (s)
S_FR = t_FR * V_TD
mu = 0.4  %Friction coefficient for braking
V1 = V_TD
V2 = 0
K_A_landing = (rho / (2 * W_S)) * (mu * CL - CD_0 - (CL^2) / (pi * AR * e))
K_T_landing = (-0.4 * T) / W_landing - mu
S_B = (1 / (2 * 9.81 * K_A_landing)) * log((K_T_landing + K_A_landing * V2^2) / (K_T_landing + K_A_landing * V1^2))

%Total Landing Distance
S_landing_total = 1.666 * (S_a + S_F + S_FR + S_B)

% Display Results
fprintf('Takeoff Ground Roll Distance (S_g): %.2f m\n', S_g);
fprintf('Takeoff Rotation Distance (S_R): %.2f m\n', S_R);
fprintf('Takeoff Transition Distance (S_TR): %.2f m\n', S_TR);
fprintf('Takeoff Climb Distance (S_CL): %.2f m\n', S_CL);
fprintf('Total Takeoff Distance (S_TO): %.2f m\n', S_TO);
fprintf('Balanced Field Length (BFL): %.2f m\n', BFL);
fprintf('Landing Approach Distance (S_a): %.2f m\n', S_a);
fprintf('Landing Flare Distance (S_F): %.2f m\n', S_F);
fprintf('Landing Free Roll Distance (S_FR): %.2f m\n', S_FR);
fprintf('Landing Braking Distance (S_B): %.2f m\n', S_B);
fprintf('Total Landing Distance: %.2f m\n', S_landing_total);

