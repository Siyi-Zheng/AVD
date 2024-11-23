clear
clc
close all
% longitudinal static stability

mac = 7.4; % aircraft mean aero chord
xcg = 39.40; % aircraft cg
claw = 5.581; % wing lift curve slope
xacw = 39.6; % aero centre of wing
kf = 1.4; % some bullshit constant
lf = 77.82; % fuselage length
wf = 6.34; % fuselage width
sw = 482; % wing area
etah = 0.9; % tailplane efficiency
clah = 4.39; % tailplane lift curve slope
AR = 8.77; % wing aspect ratio
lam = 0.25; % wing taper ratio
w_twist = 3;
b = 65; % wingspan
hh = 1.55; % vertical position of hstab
xach = 74.1; % aero centre of hstab
lh = xach - xacw; % distance between wing c/4 and hstab c/4
sweep = 26.6; % quarter chord sweep
sh = 68; % hstab area

Cmoairf = -0.131; %X-FOIL incompressible airfoil zero lift pitching moment
compressibility_factor_cruise = 1.350; %compressibility factor at cruise conditon
iw = 0;
compressibility_factor_landing = 1.024;
compressibility_factor_takeoff = 1.017;
density_cruise = 0.38;
cruise_velocity = 254;
S_w = sw; % why the fuck did you do define S_w as 440 and sw as 482 lol
Cd = 0.0161 + 0.012 + 0.003; %cd at cruise conditions of aircraft
alpha_0_w = -4.8 * pi / 180; %X-foil
alpha_0_h = 1.8 * pi / 180; % its upside-down

%%%%%%%%%%%%%%%%%%%%%%%%%%cruise
deda_cruise = (4.44 * ((1 ./ AR - 1 / (1 + AR ^ 1.7)) ...
    * ((10 - 3 * lam) / 7) * (1 - hh / b) ./ ...
    (2 * lh / b) .^ (1/3) * ...
    (cosd(sweep) .^ 0.5)) .^ 1.19) .* compressibility_factor_cruise; %downwash derivative

cmaf = kf * lf * wf ^ 2 / (mac * sw); % fuselage pitching moment

xnp = mac .* ((claw .* xacw ./ mac - cmaf + etah .* clah .* ...
    (1 - deda_cruise) .* sh ./ sw .* xach ./ mac) ...
    ./ (claw + etah .* clah .* (1 - deda_cruise) .* sh ./ sw)); % neutral pt

kn = (xnp - xcg) ./ mac; % we want this to be ~ 7 % (power off)

figure
plot(xacw, kn)

kn_on = kn - 0.02; %power on

%this must be negative to ensure longitudinal pitch stabilty
dcmcg_dalpha = (-claw .* (xacw - xcg) / mac + cmaf - etah * clah *...
    (1 - deda_cruise) * (sh / sw) * (xach / mac)) / (claw + etah * clah...
    * (1 - deda_cruise) * (sh / sw));


%%%%%%%%%%%%%%%%%%%%%%%%%%takeoff
deda_takeoff = (4.44 * ((1 / AR - 1 / (1 + AR ^ 1.7)) ...
    * ((10 - 3 * lam) / 7) * (1 - hh / b) / ...
    (2 * lh / b) ^ (1/3) * ...
    (cosd(sweep) ^ 0.5)) ^ 1.19) * compressibility_factor_takeoff; %downwash derivative

cmaf = kf * lf * wf ^ 2 / (mac * sw); % fuselage pitching moment

xnp = mac * ((claw * xacw / mac - cmaf + etah * clah * ...
    (1 - deda_takeoff) * sh / sw * xach / mac) ...
    / (claw + etah * clah * (1 - deda_takeoff) * sh / sw)); % neutral pt

kn = (xnp - xcg) / mac; % we want this to be ~7 % (power off)

kn_on = kn - 0.02; %power on

%this must be negative to ensure longitudinal pitch stabilty
dcmcg_dalpha = (-claw * (xacw - xcg) / mac + cmaf - etah * clah *...
    (1 - deda_cruise) * (sh / sw) * (xach / mac)) / (claw + etah * clah...
    * (1 - deda_cruise) * (sh / sw));

% trim analysis (clean (a/c))

CMow = (Cmoairf * (AR * cosd(sweep) ^ 2 / (AR + 2 * cosd(sweep))) -0.01...
    * w_twist) * compressibility_factor_cruise; %steady level flight
q = density_cruise * cruise_velocity ^ 2 * 0.5;
Thrust = q * S_w * Cd;

figure; % create plot of CL vs CMcg

for ih = -3:1:15
    % convert to radians
    ih = ih * pi / 180;

    temp_list_cl = []; % stores values for plotting
    temp_list_cmcg = []; % stores values for plotting

    for alpha = -5:1:20
    % i know this isnt best practice but i dont want to convert all variables to arrays
    % convert to radians
    alpha = alpha * pi / 180;

        % get coefficients
        Clw = claw * (alpha + iw - alpha_0_w);
        C_L_h = clah * ((alpha + iw - alpha_0_w) * (1 - deda_cruise)...
            + (ih - iw) - (alpha_0_h - alpha_0_w));
        C_L = Clw + etah * sh / sw * C_L_h;
        Cmcg = -Clw * ((xacw - xcg) / mac) + CMow + cmaf * alpha - etah * C_L_h...
            * (sh / S_w) * ((xach - xcg) / mac) + hh * Thrust / (q * S_w * mac);

        % plot
        temp_list_cl = [temp_list_cl, C_L];
        temp_list_cmcg = [temp_list_cmcg, Cmcg];
    end

    % plot the line
    plot(temp_list_cl, temp_list_cmcg, '-x', color="black");
    hold on;
end

% trimmed state (cruise) = ih = -9.5 (!!!), alpha = 1. In this case, CLh =
% -0.2405, CLw = 0.5650, CL = 0.5344