clear
clc

% longitudinal static stability

mac = 7.4; % aircraft mean aero chord
xcg = 37; % aircraft cg
claw = 7.74; % wing lift curve slope
xacw = 36.7; % aero centre of wing
kf = 1.4; % some bullshit constant
lf = 80.8; % fuselage length
wf = 6.34; % fuselage width
sw = 482; % wing area
etah = 0.9; % tailplane efficiency
clah = 6; % tailplane lift curve slope
AR = 8.77; % wing aspect ratio
lam = 0.25; % wing taper ratio
b = 65; % wingspan
hh = 3; % vertical position of hstab
lh = 35; % distance between wing c/4 and hstab c/4
sweep = 26.6; % quarter chord sweep
sh = 51.7; % hstab area
xach = 76; % aero centre of hstab

deda = 4.44 * ((1 / AR - 1 / (1 + AR ^ 1.7))...
    * ((10 - 3 * lam) / 7) * (1 - hh / b) /...
    (2 * lh / b) ^ (1 / 3) * ...
    (cosd(sweep) ^ 0.5)) ^ 1.19; %downwash derivative

cmaf = kf * lf * wf ^ 2 / (mac * sw); % fuselage pitching moment

xnp = mac * ((claw * xacw / mac - cmaf + etah * clah * ...
    (1 - deda) * sh / sw * xach / mac) ...
    / (claw + etah * clah * (1 - deda) * sh / sw)); % neutral pt

kn = (xnp - xcg) / mac; % we want this to be ~7% (power off)

kn_on = kn-0.02 ; %power on

%this must be negative to ensure longitudinal pitch stabilty
dcmcg_dalpha = (-claw*(xacw - xcg)/mac + cmaf - etah*clah *(1- deda)*(sh/sw)*(xach/mac))/(claw + etah*clah*(1-deda)*(sh/sw));

% trim analysis

% trim analysis