clc
clear

% Weight & Balance Estimation for Transport Aircraft

% Variable Definitions:
A = 8.77; %- Wing aspect ratio (unitless)
Ah = 5.8; %- Horizontal tailplane aspect ratio (unitless)
Av = 5.8; %- Vertical tailplane aspect ratio (unitless)
Afc = 3887.9244; %- Cargo hold floor area (ft^2), set to 0 initially, replace with actual value
Bw = 213.2; %- Wing span (ft), set to 0 initially, replace with actual value
Bh = 51; % - Horizontal tailplane span (ft)
D = 20.8; %- Maximum fuselage diameter (ft), set to 0 initially, replace with actual value
Fw = 14.7638; %- Fuselage width at horizontal tail intersection (ft), set to 0, replace as needed
HtHv = 0; % - Location of horizontal tailplane on vertical tail (0 for fuselage-mounted, 1 for T-tail)
Kbuf = 5.68; %- Factor for range; 1.02 for short ranges, 5.68 for very long ranges
Kdoor = 1; %- Door factor, 1.0 if no cargo door, 1.06 for one side cargo door, 1.12 for aft clamshell
Klav = 1.11; %- Lavatory factor, 1.11 for long range, 0.31 for short range, 3.9 for business jets
KLg = 1.12; % - Landing gear factor, 1.12 if fuselage mounted, 1.0 otherwise
Kmp = 1; %- Kneeling main gear factor, 1.126 for kneeling gear, 1.0 otherwise
Kng = 1.017; %- Nacelle factor, 1.017 for pylon mounted nacelle, 1.0 otherwise
Knp = 1; %- Nose gear factor, 1.15 for kneeling nose gear, 1.0 otherwise
Kr = 1; %      - Engine factor, 1.133 for reciprocating engines, 1.0 otherwise
Ktp = 1; %- Turboprop factor, 0.793 for turboprop, 1.0 otherwise
Kuht = 1; %- Factor for all-moving tail, 1.143 for all-moving, 1.0 otherwise
L = 246.063; %- Fuselage structural length (ft), to be assigned
La = 250; %       - Electrical routing distance; generators to avionics to cockpit (ft)
Lec = 800; %        - Engine controls routing distance; engine to cockpit - total if multi-engine (ft)
Lf = 264.4357; %- Total fuselage length (ft)
Lht = 121.4; %    - Length from wing aerodynamic center to horizontal tail aerodynamic center (ft)
Lm = 160; %         - Main landing gear length (inches)
Ln = 120; %         - Nose landing gear length (inches)
Lvt = 114.8; %    - Length from wing aerodynamic center to vertical tail aerodynamic center (ft)
Ky = 0.238 * Lht; % - Aircraft pitching radius of gyration, approx. 0.3*Lht (ft)
Kz = Lvt; % - Aircraft yaw radius of gyration, approx. Lvt (ft)
Nc = 12; %       - Number of crew, typically set as required by the design
Nen = 4; %     - Number of engines
Nf = 5; %     - Number of control functions, typically 4-7
Ngear = 2.85; % - Gear load ratio, typically 2.7-3 for commercial aircraft
Ngen = 5; %    - Number of generators, typically equal to Nen
Nl = 1.5 * Ngear; %    - Ultimate landing gear load factor, typically 1.5*Ngear
NLt = 17.4; %- Nacelle length (ft)
Nm = 1; % - Number of mechanical functions performed by controls, usually 0-2
Nmss = 4; %       - Number of main gear shock struts
Nmw = 16; %     - Number of main wheels
Nnw = 2; %     - Number of nose wheels
Np = 516; % - Total number of persons onboard (crew + passengers)
Nseat = 516; %- Number of seats of given type
Nt = 5; % - Total number of fuel tanks
Nw = 6.9; % - Nacelle width (ft)
Nz = 1.5*2.5; %- Ultimate load factor, typically 1.5 * limit load factor
Rkva = 50; % - Electrical system rating, usually 40-60 kVA for transports
Scs = 540; %    - Total control surface area (ft^2)
Scsw = 530; %   - Area of wing-mounted control surfaces (ft^2)
Se = 100; %    - Elevator area (ft^2)
Sf = 37027; %- Fuselage wetted area (ft^2)
Sht = 1039.8; %     - Horizontal tailplane area (ft^2)
Sn = 650; %     - Nacelle wetted area (ft^2)
Sw = 5188.2; %- Reference wing area (ft^2)
Svt = 433.6; %      - Vertical tailplane area (ft^2)
tc_root = 0.14; %- Wing root thickness-to-chord rEngineWenatio
tc_rootv = 0.12; %- Vertical tailplane root thickness-to-chord ratio
Vi = 59438.7; %- Integral fuel tank volume (gal)
Vp = 0; % - Self-sealing tank volume (gal)
Vpr = 150000; %     - Volume of pressurized sections (ft^3)
Vs = 231.6; %- Landing stall speed (ft/s)
Vt = Vi +Vp; %- Total volume of fuel tanks (gal)
WAPU = 500; %- Uninstalled APU weight (lb)
Wc = 0; % - Maximum cargo weight (lb)
Wdg = 859928.486; % - Design gross weight (lb)
Wen = 13552; %- Engine weight (lb)
Wenc = 19841.6; %- Weight of engine and contents (lb)
Wl = 859928.486; %- Landing design gross weight (lb)
Wseat = 13; %- Weight of single seat (lb)
Wuav = 1100; % - Uninstalled avionics weight, typically 800-1400 lb
lambda = 0.25; %- Wing taper ratio
Lambda = 0.464; %- Wing quarter chord sweep angle (radians)
Lambda_ht = 0.26; %- Horizontal tailplane quarter chord sweep (radians)
Lambda_vt = 0.26; %- Vertical tailplane quarter chord sweep (radians)
g = 32.18504;
Iy = Wdg * Ky ^ 2 / g; %- Pitching moment of inertia (lb*ft^2), approximated as Wo*Ky^2/g
Kws = 0.75 * ((1 + 2 * lambda) / (1 + lambda)) * Bw * tan(Lambda) / L; % - Wing sweep factor, calculated as 0.75 * [(1 + 2*lambda)/(1 + lambda)] * Bw * tan(Lambda) / L

% % Default values for unspecified variables (set to zero, modify as needed):
% A = 0; % Wing aspect ratio
% Ah = 0; % Horizontal tailplane aspect ratio
% Av = 0; % Vertical tailplane aspect ratio
% Bw = 0; % Wing span in ft
% D = 0; % Fuselage diameter in ft
% % ... Continue initializing all variables used in equations as above

% Begin Calculations
% 1. Aircraft Wing Weight
% Wdg = 0; % Design gross weight, replace with actual value
% Nz = 0; % Ultimate load factor, typically 1.5 * limit load factor
% Sw = 0; % Reference wing area, replace with actual value
% A = 0; % Aspect ratio of the wing
% lambda = 0; % Wing taper ratio
% Scsw = 0; % Area of wing-mounted control surfaces
% Lambda = 0; % Wing quarter chord sweep (radians)
% tc_root = 0; % Wing root thickness-to-chord ratio


% Equation for Wing Weight (Ww)
CGw = 39.6; % Wing CG location (m)
Zw = -1.83;
Ww = (0.78 * 0.0051 * (Wdg * Nz) ^ 0.557 * Sw ^ 0.649 * A ^ 0.5 * (1 + lambda) ^ 0.1 * ...
    Scsw ^ 0.1 )/ (cos(Lambda) * (tc_root) ^ 0.4); % Wing Weight
% Comment: Replace variables with actual values as needed.

% engine weight
CGeng = CGw - 0.7;
Zeng = -4.5;
W_engines = Nen * Wen;

% 2. Horizontal Tailplane Weight (Wht)
% Kuht = 1; % Factor for all-moving tail, typically 1.143 for all-moving, 1.0 otherwise
% Sht = 0; % Horizontal tailplane area in ft^2
% Ky = 0; % Aspect ratio of the horizontal tail
% Fw = 0; % Fuselage width at horizontal tail intersection
% Bh = 0; % Horizontal tailplane span in ft
% Lht = 0; % Length from wing aerodynamic center to horizontal tail aerodynamic center
% Lambda_ht = 0; % Horizontal tail quarter-chord sweep (radians)
CGht = 72.1;
Zht = 1.55;
Wht = 0.75 * 0.0379 * Kuht * Wdg ^ 0.639 * Nz ^ 0.1 * Sht ^ 0.75 * Ky ^ 0.704 * ...
    Ah ^ 0.166 * (1 + Se / Sht) ^ 0.1 / ((1 + Fw / Bh) ^ 0.25 * Lht * cos(Lambda_ht));
% Comment: Replace with actual values for Wht calculation.

% 3. Vertical Tailplane Weight (Wvt)
% HtHv = 0; % Location factor, 0.0 for fuselage mounted, 1.0 for T-tail
% Svt = 0; % Vertical tailplane area in ft^2
% Av = 0; % Vertical tailplane aspect ratio
% Lvt = 0; % Length from wing aerodynamic center to vertical tail aerodynamic center
% Lambda_vt = 0; % Vertical tail quarter-chord sweep (radians)
% tc_rootv = 0; % Vertical tail root thickness-to-chord ratio
CGvt = 71.1;
Zvt = 4.64;
Wvt = 0.75 * 0.0026 * (1 + HtHv) ^ 0.225 * Wdg ^ 0.556 * Nz ^ 0.536 * Svt ^ 0.5 * ...
    Ky ^ 0.875 * Av ^ 0.35 * Lvt ^ 0.5 * cos(Lambda_vt) * (tc_rootv) ^ 0.5;
% Comment: Replace variables for accurate vertical tail weight estimation.

% Continue defining
% Continue defining remaining equations based on the document

% 4. Fuselage Weight (Wfus)
% Sf = 0; % Fuselage wetted area (ft^2)
% L_D = 0; % Lift-to-drag ratio, typically 10-20 for transport aircraft
CGfus = L * 0.435 / 3.281;
Zfus = 0;
Wfus = 0.85 * 0.3280 * Kdoor * KLg * (Wdg * Nz) ^ 0.5 * L ^ 0.25 * Sf ^ 0.302 * ...
    (1 + Kws) ^ 0.04 * (L/D) ^ 0.1;
% Comment: Replace variables with actual values for fuselage weight.

% 5. Main Landing Gear Weight (Wmlg)
% Lm = 0; % Main landing gear length (inches)
% Nmss = 0; % Number of main gear shock struts
% Nmw = 0; % Number of main wheels
% Vs = 0; % Landing stall speed (ft/s)
CGmlg = 43;
Zmlg = -6.33;
Wmlg = 0.88 * 0.0106 * Kmp * Wl ^ 0.888 * Nl ^ 0.25* Lm ^ 0.4 *Nmw^0.321 * Vs ^ 0.1 / (Nmss ^ 0.5);
% Comment: Replace values for main landing gear weight.

% 6. Nose Landing Gear Weight (Wnlg)
% Ln = 0; % Nose landing gear length (inches)
% Nnw = 0; % Number of nose wheels
CGnlg = 7;
Znlg = -6.33;
Wnlg = 0.88 * 0.032 * Knp * Wl ^ 0.646 * Nl ^ 0.2 * Ln ^ 0.5 * Nnw ^ 0.45;
% Comment: Replace variables with actual values for nose landing gear weight.

% 7. Nacelle Weight (Winl)
% NLt = 0; % Nacelle length (ft)
% Nw = 0; % Nacelle width (ft)
% Wenc = 0; % Weight of engine and contents (lb)
% Sn = 0; % Nacelle wetted area (ft^2)
CGinl = 38.12;
Zinl = Zeng;
Winl = 0.85 * 0.6724 * Kng * NLt ^ 0.1 * Nw ^ 0.294 * Nz ^ 0.119 * Wenc ^ 0.611 *Nen^0.984* Sn ^ 0.224;
% Comment: Replace values for nacelle weight.

% 8. Engine Controls Weight (Wec)
% Lec = 0; % Engine control routing distance (ft)
CGec = CGeng;
Zec = Zeng;
Wec = 5 * Nen + 0.8 * Lec;
% Comment: Weight of engine controls; set Nen and Lec as needed.

% 9. Engine Starter Weight (Wes)
% Wen = 0; % Engine weight (lb)
CGes = CGeng;
Zes = Zeng;
Wes = 49.19 * (Nen * Wen / 1000) ^ 0.541;
% Comment: Engine starter weight; specify engine weight and number of engines.

% 10. Fuel System Weight (Wfs)
% Vt = 0; % Total volume of fuel tanks (gal)
% Nt = 0; % Number of fuel tanks
% Vp = 0; % Self-sealing tank volume (gal)
% Vi = 0; % Integral fuel tank volume (gal)
CGfs = 39.5;
Zfs = -2.01;
Wfs = 2.405 * Vt ^ 0.606 * Nt ^ 0.5 * (1 + Vp / Vt) / (1 + Vi / Vt);
% Comment: Fuel system weight; define fuel volumes and tank counts.

% 11. Flight Controls Weight (Wfc)
% Scs = 0; % Total control surface area (ft^2)
% Iy = 0; % Pitching moment of inertia (lb*ft^2), usually ≈ Wdg * Ky^2 / g
% Nm = 0; % Number of mechanical functions performed by controls
CGfc = 50;
Zfc = -1;
Wfc = 145.9 * Nf ^ 0.554 * Scs ^ 0.2 * (Iy * 1e-6) ^ 0.07 / (1 + Nm / Nf);
% Comment: Flight controls weight; specify control surface area and inertia.

% 12. Installed APU Weight (WAPUinst)
% WAPU = 0; % Uninstalled APU weight (lb)
CGapu = 75;
Zapu = 1.34;
WAPUinst = 2.2 * WAPU;
% Comment: Installed APU weight; specify WAPU.

% 13. Instruments Weight (Winstr)
% Kr = 1; % Engine type factor, 1.133 for reciprocating engines, 1.0 otherwise
% Ktp = 1; % Turboprop factor, 0.793 for turboprop, 1.0 otherwise
CGinstr = 4;
Zinstr = 0;
Winstr = 4.509 * Kr * Ktp * Nc ^ 0.541 * Nen * (Lf + Bw) ^ 0.5;
% Comment: Instrument weight; define factors, number of crew, and dimensions.

% 14. Hydraulic System Weight (Whydr)
CGhydr = 50;
Zhydr = -1;
Whydr = 0.2673 * Nf * (Lf + Bw) ^ 0.937;
% Comment: Hydraulic system weight; replace Nf, Lf, Bw as needed.

% 15. Electrical System Weight (Wel)
% Rkva = 0; % System electrical rating (kVA)
CGel = 45;
Zel = -1;
Wel = 7.291 * Rkva ^ 0.782 * La ^ 0.346 * Ngen ^ 0.1;
% Comment: Electrical system weight; set electrical rating and routing distance.

% 16. Avionics Weight (Wav)
% Wuav = 0; % Uninstalled avionics weight, typically 800-1400 lb
CGav = 5;
Zav = 0;
Wav = 1.73 * Wuav ^ 0.983;
% Comment: Avionics weight; specify Wuav.

% 17. Furnishings Weight (Wfurn)
% Nseat = 0; % Number of seats
% Wseat = 0; % Weight of single seat, e.g., 60 lb for flight deck, 32 lb for passenger
% Klav = 0; % Lavatory factor
% Kbuf = 0; % Buffer factor for furnishings
CGfurn = CGfus;
Zfurn = 1;
Wfurn = 0.0577 * Nc ^ 0.1 * Wc ^ 0.393 * Sf ^ 0.75 + Nseat * Wseat + Klav * Np ^ 1.33 + Kbuf * Np ^ 1.12;
% Comment: Furnishings weight; define seating and comfort parameters.

% 18. Air-Conditioning Weight (Wac)
% Vpr = 0; % Volume of pressurized sections (ft^3)
CGac = CGfus;
Zac = 2.5;
Wac = 62.36 * Np ^ 0.25 * (Vpr * 1e-3) ^ 0.604 * Wuav ^ 0.1;
% Comment: Air-conditioning weight; specify volumes and number of persons.

% 19. Anti-Icing System Weight (Wai)
CGai = CGw - 1;
Zai = Zw;
Wai = 0.002 * Wdg;
% Comment: Anti-icing system weight; set Wdg.

% 20. Handling Gear Weight (Whg)
CGwhg = 0;
Whg_civilian = 0; % Civilian handling gear weight
% Afc = 0; % Cargo hold floor area (ft^2)
CGmil = 35;
Whg_military = 0; % Military handling gear weight

% Final output
% Display calculated weights for verification:
disp(['Wing Weight (Ww): ', num2str(Ww*0.453592)]);
disp(['Horizontal Tailplane Weight (Wht): ', num2str(Wht*0.453592)]);
disp(['Vertical Tailplane Weight (Wvt): ', num2str(Wvt*0.453592)]);
disp(['Fuselage Weight (Wfus): ', num2str(Wfus*0.453592)]);
disp(['Main Landing Gear Weight (Wmlg): ', num2str(Wmlg*0.453592)]);
disp(['Nose Landing Gear Weight (Wnlg): ', num2str(Wnlg*0.453592)]);
disp(['Nacelle Weight (Winl): ', num2str(Winl*0.453592)]);
disp(['Engine Controls Weight (Wec): ', num2str(Wec*0.453592)]);
disp(['Engine Starter Weight (Wes): ', num2str(Wes*0.453592)]);
disp(['Fuel System Weight (Wfs): ', num2str(Wfs*0.453592)]);
disp(['Flight Controls Weight (Wfc): ', num2str(Wfc*0.453592)]);
disp(['Installed APU Weight (WAPUinst): ', num2str(WAPUinst*0.453592)]);
disp(['Instruments Weight (Winstr): ', num2str(Winstr*0.453592)]);
disp(['Hydraulic System Weight (Whydr): ', num2str(Whydr*0.453592)]);
disp(['Electrical System Weight (Wel): ', num2str(Wel*0.453592)]);
disp(['Avionics Weight (Wav): ', num2str(Wav*0.453592)]);
disp(['Furnishings Weight (Wfurn): ', num2str(Wfurn*0.453592)]);
disp(['Air-Conditioning Weight (Wac): ', num2str(Wac*0.453592)]);
disp(['Anti-Icing System Weight (Wai): ', num2str(Wai*0.453592)]);
disp(['Civilian Handling Gear Weight (Whg_civilian): ', num2str(Whg_civilian*0.453592)]);
disp(['Military Handling Gear Weight (Whg_military): ', num2str(Whg_military*0.453592)]);
disp(['Engine Weight (W_engines): ', num2str(W_engines*0.453592)]);

% get total weight
Wtotal = Ww + Wht + Wvt + Wfus + Wmlg + Wnlg + Winl + Wec + Wes + Wfs + ...
    Wfc + WAPUinst + Winstr + Whydr + Wel + Wav + Wfurn + Wac + Wai + ...
    Whg_civilian + Whg_military + W_engines;

Wtotal_tons = Wtotal / 2204;
disp(Wtotal_tons)

disp(['Total Weight (empty): ', num2str(Wtotal_tons), ' tons']);

% get cg
CGtotal = (Ww * CGw + Wht * CGht + Wvt * CGvt + Wfus * CGfus + Wmlg * CGmlg + ...
    Wnlg * CGnlg + Winl * CGinl + Wec * CGec + Wes * CGes + Wfs * CGfs + ...
    Wfc * CGfc + WAPUinst * CGapu + Winstr * CGinstr + Whydr * CGhydr + ...
    Wel * CGel + Wav * CGav + Wfurn * CGfurn + Wac * CGac + Wai * CGai + ...
    Whg_civilian * CGwhg + Whg_military * CGmil + W_engines * CGeng) / Wtotal;

disp(['Total x_cg (empty): ', num2str(CGtotal), ' m']);

% get z center of gravity

ZCGtotal = (Ww * Zw + Wht * Zht + Wvt * Zvt + Wfus * Zfus + Wmlg * Zmlg + ...
    Wnlg * Znlg + Winl * Zinl + Wec * Zec + Wes * Zes + Wfs * Zfs + ...
    Wfc * Zfc + WAPUinst * Zapu + Winstr * Zinstr + Whydr * Zhydr + ...
    Wel * Zel + Wav * Zav + Wfurn * Zfurn + Wac * Zac + Wai * Zai) / Wtotal;

disp(['Total z_cg (empty): ', num2str(ZCGtotal), ' m']);

% add fuel, pax, luggage
Wfuel = 170; % tons
Wpax = 38.7; % tons
Wluggage = 9.42; % tons

Wtotal_full = Wtotal_tons + Wfuel + Wpax + Wluggage;
disp(['Total Weight (full): ', num2str(Wtotal_full), ' tons']);


% get cg
CGfuel = CGfs; % m
CGpax = 40; % m
CGluggage = 40; % m
Zfuel = Zfs;
Zpax = 0.5;
Zluggage = -0.5;

CGtotal_full = (Wtotal_tons * CGtotal + Wfuel * CGfuel + Wpax * CGpax + Wluggage * CGluggage) / Wtotal_full;
disp(['Total CG (full): ', num2str(CGtotal_full), ' m']);


zCGtotal_full = (Wtotal_tons * ZCGtotal + Wfuel * Zfuel + Wpax * Zpax + Wluggage * Zluggage) / Wtotal_full;
disp(['Total z_cg (full): ', num2str(zCGtotal_full), ' m']);

%plotting the CG envelope

W_cruise_start= Wtotal_tons + 0.97*0.985*Wfuel + Wpax + Wluggage;
CG_cruise_start= (Wtotal_tons * CGtotal + 0.97*0.985*Wfuel * CGfuel + Wpax * CGpax + Wluggage * CGluggage) / W_cruise_start;
disp(['Total CG (cruise start): ', num2str(CG_cruise_start), ' m']);

W_cruise_end= Wtotal_tons + 0.97*0.985*0.6225*Wfuel + Wpax + Wluggage;
CG_cruise_end= (Wtotal_tons * CGtotal + 0.97*0.985*0.625*Wfuel * CGfuel + Wpax * CGpax + Wluggage * CGluggage) / W_cruise_end;
disp(['Total CG (cruise end): ', num2str(CG_cruise_end), ' m']);

Wtotal_nofuel= Wtotal_tons + Wpax + Wluggage;
disp(['Total Weight (no fuel): ', num2str(Wtotal_nofuel), ' tons']);
CGtotal_nofuel = (Wtotal_tons * CGtotal + Wpax * CGpax + Wluggage * CGluggage) / Wtotal_nofuel;
disp(['Total CG (no fuel): ', num2str(CGtotal_nofuel), ' m']);





figure
scatter(CG_cruise_start, W_cruise_start, 'green')
hold on
scatter(CG_cruise_end, W_cruise_end, 'blue')
scatter(CGtotal_nofuel, Wtotal_nofuel, 'red')
scatter(CGtotal_full, Wtotal_full, 'magenta')


% get Iyy using sum(Iyy) = sum(m * (x - x_cg)^2)
Iyy = Ww * (CGw - CGtotal)^2 + Wht * (CGht - CGtotal)^2 + Wvt * (CGvt - CGtotal)^2 + ...
    Wfus * (CGfus - CGtotal)^2 + Wmlg * (CGmlg - CGtotal)^2 + Wnlg * (CGnlg - CGtotal)^2 + ...
    Winl * (CGinl - CGtotal)^2 + Wec * (CGec - CGtotal)^2 + Wes * (CGes - CGtotal)^2 + ...
    Wfs * (CGfs - CGtotal)^2 + Wfc * (CGfc - CGtotal)^2 + WAPUinst * (CGapu - CGtotal)^2 + ...
    Winstr * (CGinstr - CGtotal)^2 + Whydr * (CGhydr - CGtotal)^2 + Wel * (CGel - CGtotal)^2 + ...
    Wav * (CGav - CGtotal)^2 + Wfurn * (CGfurn - CGtotal)^2 + Wac * (CGac - CGtotal)^2 + ...
    Wai * (CGai - CGtotal)^2 + Whg_civilian * (CGwhg - CGtotal)^2 + Whg_military * (CGmil - CGtotal)^2 + ...
    W_engines * (CGeng - CGtotal)^2 + Wfuel * (CGfuel - CGtotal_full)^2 + ...
    Wpax * (CGpax - CGtotal_full)^2 + Wluggage * (CGluggage - CGtotal_full)^2;

% convert from lb*m^2 to kg*m^2
Iyy = Iyy * 0.453592;

disp(['Iyy: ', num2str(Iyy), ' kg*m^2']);