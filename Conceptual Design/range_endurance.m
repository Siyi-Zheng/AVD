clc
clear

% we have to keep mwf same as before
MWF = [0.9700    0.9850    0.6225    0.9900    0.9850    0.9873    0.9844    0.9900    0.9950];

% get (L/D)_max from updated aero calcs
L_D_max = 19.0692; % i think this was right???
L_D_cruise = 0.866 * L_D_max;
L_D_loiter = L_D_max;

% sfc values
sfc_cruise = 13.825 * 1e-6 * 9.81; % 1/s
sfc_loiter = 11.1 * 1e-6 * 9.81; % 1/s

% velocity values
V_cruise1 = 245; % m/s
V_cruise2 = 241.5; % m/s
V_loiter = 120; % m/s

% range calculation (cruise 1)
% get relevant mass weight fractions
mwf_cruise1 = MWF(3);
mwf_cruise2 = MWF(6);
mwf_loiter = MWF(7);
% breguet range equation
range_cruise1 = V_cruise1 * (1/sfc_cruise) * L_D_cruise * log(1/mwf_cruise1);
range_cruise2 = V_cruise2 * (1/sfc_cruise) * L_D_cruise * log(1/mwf_cruise2);
endurance_loiter = (1/sfc_loiter) * L_D_loiter * log(1/mwf_loiter);

% do the payload graph thing

% max fuel, zero payload
W_max = 353385;
MTOW = 390000;
W_empty = 142971; % kg
W_fuel_cruise = W_max * prod(MWF(1:2)) * (1 - MWF(3));
W_post_cruise = W_empty / prod(MWF(4:end));
W_pre_cruise = W_post_cruise + W_fuel_cruise;
range_fuelonly = max(range_cruise1, V_cruise1 * (1/sfc_cruise) * L_D_cruise * ...
    log(W_pre_cruise/W_post_cruise)) / 1852;

% max fuel, MTOW
W_payload = W_max - W_empty / prod(MWF);
W_post_cruise = (W_empty + W_payload) / prod(MWF(4:end));
W_pre_cruise = W_post_cruise + W_fuel_cruise;
range_maxfuel = max(range_cruise1, V_cruise1 * (1/sfc_cruise) * L_D_cruise * ...
    log(W_pre_cruise/W_post_cruise)) / 1852;

% max payload, MTOW
W_payload2 = MTOW - W_empty / prod(MWF);
W_post_cruise2 = (W_empty + W_payload2) / prod(MWF(4:end));
W_pre_cruise2 = W_post_cruise2 + W_fuel_cruise;
range_maxpayload = V_cruise1 * (1/sfc_cruise) * L_D_cruise * ...
    log(W_pre_cruise2/W_post_cruise2) / 1852;

% range plot
figure;
plot([0 range_maxpayload], [W_payload2 W_payload2] / 1000, "k")
hold on
plot([range_maxfuel range_fuelonly], [W_payload 0]/1000, "k")
plot([range_maxpayload range_maxfuel], [W_payload2 W_payload]/1000, "k")
xline(7500, "--")
xlabel("Cruise Range (nmi)")
ylabel("Payload mass (tonnes)")
legend("Aircraft Cruise Range Envelope", "", "Cruise Range Requirement")
grid on