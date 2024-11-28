clc
clear
close all

% Weight & Balance Estimation for Transport Aircraft

% Variable Definitions:
mac=7.4;
A = 8.77; %- Wing aspect ratio (unitless)
Ah = 5.8; %- Horizontal tailplane aspect ratio (unitless)
Av = 2*1.6; %- Vertical tailplane aspect ratio (unitless)
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
La = 404; %       - Electrical routing distance; generators to avionics to cockpit (ft)
Lec = 266; %        - Engine controls routing distance; engine to cockpit - total if multi-engine (ft)
Lf = 264.4357; %- Total fuselage length (ft)
Lht = 106.6; %    - Length from wing aerodynamic center to horizontal tail aerodynamic center (ft)
Lm = 160; %         - Main landing gear length (inches)
Ln = 120; %         - Nose landing gear length (inches)
Lvt = 103.6; %    - Length from wing aerodynamic center to vertical tail aerodynamic center (ft)
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
tc_root = 0.14; %- Wing root thickness-to-chord ratio
tc_rootv = 0.12; %- Vertical tailplane root thickness-to-chord ratio
Vi = 59438.7; %- Integral fuel tank volume (gal)
Vp = 0; % - Self-sealing tank volume (gal)
Vpr = 150000; %     - Volume of pressurized sections (ft^3)
Vs = 231.6; %- Landing stall speed (ft/s)
Vt = Vi +Vp; %- Total volume of fuel tanks (gal)
WAPU = 500; %- Uninstalled APU weight (lb)
Wc = 0; % - Maximum cargo weight (lb)
Wdg = 779200; % - Design gross weight (lb)
Wen = 13552; %- Engine weight (lb)
Wenc = 19841.6; %- Weight of engine and contents (lb)
Wl = Wdg * 0.85; %- Landing design gross weight (lb)
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
CGw = 35; % Wing CG location (m)
Zw = -1.83;
Ww = 1.05 * (0.78 * 0.0051 * (Wdg * Nz) ^ 0.557 * Sw ^ 0.649 * A ^ 0.5 * (1 + lambda) ^ 0.1 * ...
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
Wht = 1.05 * 0.75 * 0.0379 * Kuht * Wdg ^ 0.639 * Nz ^ 0.1 * Sht ^ 0.75 * Ky ^ 0.704 * ...
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
Wvt = 1.05 * 0.75 * 0.0026 * (1 + HtHv) ^ 0.225 * Wdg ^ 0.556 * Nz ^ 0.536 * Svt ^ 0.5 * ...
    Kz ^ 0.875 * Av ^ 0.35 / (Lvt ^ 0.5 * cos(Lambda_vt) * (tc_rootv) ^ 0.5);
% Comment: Replace variables for accurate vertical tail weight estimation.

% Continue defining
% Continue defining remaining equations based on the document

% 4. Fuselage Weight (Wfus)
% Sf = 0; % Fuselage wetted area (ft^2)
% L_D = 0; % Lift-to-drag ratio, typically 10-20 for transport aircraft
CGfus = L * 0.435 / 3.281;
Zfus = 0;
Wfus = 1.05 * 0.85 * 0.3280 * Kdoor * KLg * (Wdg * Nz) ^ 0.5 * L ^ 0.25 * Sf ^ 0.302 * ...
    (1 + Kws) ^ 0.04 * (L/D) ^ 0.1;
% Comment: Replace variables with actual values for fuselage weight.

% 5. Main Landing Gear Weight (Wmlg)
% Lm = 0; % Main landing gear length (inches)
% Nmss = 0; % Number of main gear shock struts
% Nmw = 0; % Number of main wheels
% Vs = 0; % Landing stall speed (ft/s)
CGmlg = 43;
Zmlg = -6.33;
Wmlg = 1.05 * 0.88 * 0.0106 * Kmp * Wl ^ 0.888 * Nl ^ 0.25 * Lm ^ 0.4 *Nmw^0.321 * Vs ^ 0.1 / (Nmss ^ 0.5);
% Comment: Replace values for main landing gear weight.

% 6. Nose Landing Gear Weight (Wnlg)
% Ln = 0; % Nose landing gear length (inches)
% Nnw = 0; % Number of nose wheels
CGnlg = 7;
Znlg = -6.33;
Wnlg = 1.05 * 0.88 * 0.032 * Knp * Wl ^ 0.646 * Nl ^ 0.2 * Ln ^ 0.5 * Nnw ^ 0.45;
% Comment: Replace variables with actual values for nose landing gear weight.

% 7. Nacelle Weight (Winl)
% NLt = 0; % Nacelle length (ft)
% Nw = 0; % Nacelle width (ft)
% Wenc = 0; % Weight of engine and contents (lb)
% Sn = 0; % Nacelle wetted area (ft^2)
CGinl = CGw - 1.4;
Zinl = Zeng;
Winl = 1.05 * 0.85 * 0.6724 * Kng * NLt ^ 0.1 * Nw ^ 0.294 * Nz ^ 0.119 * Wenc ^ 0.611 *Nen^0.984* Sn ^ 0.224;
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
CGfs = CGeng + 0.1;
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

% add thrust reverser effects
W_engines = W_engines + 0.3 * Winl;

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
Wfuel = 176; % tons
Wpax = 38.7; % tons
Wluggage = 9.42; % tons

Wtotal_full = Wtotal_tons + Wfuel + Wpax + Wluggage;
disp(['Total Weight (full): ', num2str(Wtotal_full), ' tons']);


% get cg
CGfuel = CGfs; % m
CGpax = CGfus; % m
CGluggage = CGfus; % m
Zfuel = Zfs;
Zpax = 0.5;
Zluggage = -0.5;


%==================== FUEL TANK DISTRIBUTION ==============================
CG_iwtank = CGw - 39.6 + 36.553; % inner wing tank xcg (m)
CG_ftank = CGw - 39.6 + 33.39 + 5; % fuselage tank xcg (m)
CG_owtank = CGw - 39.6 + 42.01; % outer wing tank xcg (m)
CG_trimtank = 73; % trim tank xcg (m)

W_iwtank = 55.5; % inner wing tank fuel mass (tons)
W_ftank = 44.2; % fuselage tank fuel mass (tons)
W_owtank = 76.4; % outer wing tank fuel mass (tons)
W_trimtank = 0; % trim tank fuel mass (tons)

CG_fore = (CG_iwtank * W_iwtank + CG_ftank * W_ftank) / (W_iwtank + W_ftank);
W_fore = (W_iwtank + W_ftank);

CGfuel = (CG_iwtank * W_iwtank + CG_ftank * W_ftank + CG_owtank *...
    W_owtank + CG_trimtank * W_trimtank) / (W_iwtank + W_ftank...
    + W_owtank + W_trimtank);

CG_fuelled_front = (Wtotal_tons * CGtotal + W_fore * CG_fore + Wpax...
    * CGpax + Wluggage * CGluggage) / (Wtotal_tons + W_fore + Wpax + Wluggage);

CG_fuelled_aft = (Wtotal_tons * CGtotal + W_owtank * CG_owtank + Wpax...
    * CGpax + Wluggage * CGluggage) / (Wtotal_tons + W_owtank + Wpax + Wluggage);

CGtotal_full = (Wtotal_tons * CGtotal + (Wfuel+W_trimtank) .* CGfuel + Wpax * CGpax...
    + Wluggage * CGluggage) ./ (Wtotal_full + W_trimtank);
disp(['Total CG (full): ', num2str(CGtotal_full), ' m']);

CGempty_trimtank = (Wtotal_tons * CGtotal + W_trimtank * CG_ftank)...
    / (Wtotal_tons + W_trimtank);

zCGtotal_full = (Wtotal_tons * ZCGtotal + Wfuel * Zfuel + Wpax * Zpax...
    + Wluggage * Zluggage) / (Wtotal_full + W_trimtank);
disp(['Total z_cg (full): ', num2str(zCGtotal_full), ' m']);

%plotting the CG envelope

%%cruise start
W_cruise_start= Wtotal_tons + 0.97*0.985*Wfuel + Wpax + Wluggage;
CG_cruise_start= (Wtotal_tons * CGtotal + 0.97*0.985*Wfuel * CGfuel + Wpax...
    * CGpax + Wluggage * CGluggage) / W_cruise_start;
disp(['Total CG (cruise start): ', num2str(CG_cruise_start), ' m']); %operational region

%%with passengers
W_owtank_remaining=  W_owtank - (Wfuel-( 0.97*0.985* Wfuel));
W_cruise_start_fore= Wtotal_tons + W_owtank_remaining + Wpax + Wluggage + W_trimtank + W_fore;
CG_cruise_start_fore= (Wtotal_tons * CGtotal + (W_owtank_remaining * CG_owtank) + W_iwtank * CG_iwtank +(W_ftank * CG_ftank) + Wpax * CGpax + Wluggage * CGluggage) / W_cruise_start_fore;
disp(['Total CG (cruise start fore): ', num2str(CG_cruise_start_fore), ' m']); %assuming asymmetric fuel burn

W_iwtank_remaining= W_iwtank - (Wfuel-(0.97*0.985* Wfuel));
W_cruise_start_aft = Wtotal_tons + Wpax + Wluggage + W_trimtank + W_iwtank_remaining + W_owtank + W_ftank;
CG_cruise_start_aft= (Wtotal_tons * CGtotal + W_iwtank_remaining*CG_iwtank + W_owtank * CG_owtank + W_ftank * CG_ftank + Wpax * CGpax + Wluggage * CGluggage) / W_cruise_start_aft;
disp(['Total CG (cruise start aft): ', num2str(CG_cruise_start_aft), ' m']); %assuming asymmetric fuel burn

%%without passengers
W_owtank_remaining=  W_owtank - (Wfuel-( 0.97*0.985* Wfuel));
W_cruise_start_fore_nopax= Wtotal_tons + W_owtank_remaining + W_trimtank + W_fore;
CG_cruise_start_fore_nopax= (Wtotal_tons * CGtotal + (W_owtank_remaining * CG_owtank) + W_iwtank * CG_iwtank +(W_ftank * CG_ftank) )/ W_cruise_start_fore_nopax;
disp(['Total CG (cruise start fore_nopax): ', num2str(CG_cruise_start_fore_nopax), ' m']); %assuming asymmetric fuel burn

W_iwtank_remaining= W_iwtank - (Wfuel-(0.97*0.985* Wfuel));
W_cruise_start_aft_nopax = Wtotal_tons + W_trimtank + W_iwtank_remaining + W_owtank + W_ftank;
CG_cruise_start_aft_nopax= (Wtotal_tons * CGtotal + W_iwtank_remaining*CG_iwtank + W_owtank * CG_owtank + W_ftank * CG_ftank) / W_cruise_start_aft_nopax;
disp(['Total CG (cruise start aft_nopax): ', num2str(CG_cruise_start_aft_nopax), ' m']); %assuming asymmetric fuel burn


%%End of cruise
W_cruise_end= Wtotal_tons + 0.97*0.985*0.6225*Wfuel + Wpax + Wluggage;
CG_cruise_end= (Wtotal_tons * CGtotal + 0.97*0.985*0.6225*Wfuel * CGfuel...
    + Wpax * CGpax + Wluggage * CGluggage) / W_cruise_end;
disp(['Total CG (cruise end): ', num2str(CG_cruise_end), ' m']);

%%with passengers
W_owtank_leftend= W_owtank - (Wfuel-(0.97*0.985*0.6225*Wfuel)); %all outerwing used up
W_cruise_end_fore= Wtotal_tons + W_owtank_leftend + Wpax + Wluggage+ W_trimtank +W_fore ;
CG_cruise_end_fore= (Wtotal_tons * CGtotal + W_owtank_leftend*CG_owtank + W_iwtank*CG_iwtank + W_ftank*CG_ftank + Wpax * CGpax + Wluggage * CGluggage) / W_cruise_end_fore;
disp(['Total CG (cruise end fore): ', num2str(CG_cruise_end_fore), ' m']);

W_iwtank_leftend= 0;
W_ftank_leftend = abs(W_iwtank - (Wfuel-(0.97*0.985*0.6225* Wfuel))) ;
W_cruise_end_aft = Wtotal_tons + Wpax + Wluggage + W_trimtank + W_iwtank_leftend + W_owtank + W_ftank_leftend;
CG_cruise_end_aft= (Wtotal_tons * CGtotal + W_iwtank_leftend*CG_iwtank + W_owtank * CG_owtank + W_ftank_leftend * CG_ftank + Wpax * CGpax + Wluggage * CGluggage) / W_cruise_end_aft;
disp(['Total CG (cruise end aft): ', num2str(CG_cruise_end_aft), ' m']); %assuming asymmetric fuel burn

%%without passengers
W_owtank_leftend_nopax= W_owtank - (Wfuel-(0.97*0.985*0.6225*Wfuel)); %all outerwing used up
W_cruise_end_fore_nopax= Wtotal_tons + W_owtank_leftend_nopax + W_trimtank +W_fore ;
CG_cruise_end_fore_nopax= (Wtotal_tons * CGtotal + W_owtank_leftend*CG_owtank + W_iwtank*CG_iwtank + W_ftank*CG_ftank) / W_cruise_end_fore_nopax;
disp(['Total CG (cruise end fore no pax): ', num2str(CG_cruise_end_fore_nopax), ' m']);

W_iwtank_leftend_nopax= 0;
W_ftank_leftend_nopax= abs(W_iwtank - (Wfuel-(0.97*0.985*0.6225* Wfuel))) ;
W_cruise_end_aft_nopax = Wtotal_tons + W_trimtank + W_iwtank_leftend_nopax + W_owtank + W_ftank_leftend_nopax;
CG_cruise_end_aft_nopax= (Wtotal_tons * CGtotal + W_iwtank_leftend_nopax*CG_iwtank + W_owtank * CG_owtank + W_ftank_leftend_nopax * CG_ftank) / W_cruise_end_aft_nopax;
disp(['Total CG (cruise end aft no pax): ', num2str(CG_cruise_end_aft_nopax), ' m']); %assuming asymmetric fuel burn

Wtotal_nofuel= Wtotal_tons + Wpax + Wluggage; 
disp(['Total Weight (no fuel): ', num2str(Wtotal_nofuel), ' tons']);
CGtotal_nofuel = (Wtotal_tons * CGtotal + Wpax * CGpax + Wluggage * CGluggage) / Wtotal_nofuel;
disp(['Total CG (no fuel): ', num2str(CGtotal_nofuel), ' m']);

%%LANDING 

W_landing = Wtotal_tons + (0.9700* 0.9850 *0.6225* 0.9900*0.9850* 0.9873*0.9844* 0.9900* 0.9950)*Wfuel + Wpax + Wluggage ;
CG_landing_design= (Wtotal_tons * CGtotal + 0.97*0.985*0.6225* 0.9900*0.9850* 0.9873*0.9844* 0.9900* 0.9950*Wfuel * CGfuel + Wpax...
    * CGpax + Wluggage * CGluggage) / W_landing;
disp(['CG landing: ', num2str(CG_landing_design), ' m']);

%%with passengers
W_owtank_remaining_landing=0;
W_iwtank_remaining_landing= abs(W_owtank - (Wfuel - (0.9700* 0.9850 *0.6225* 0.9900*0.9850* 0.9873*0.9844* 0.9900* 0.9950)*Wfuel)) ;
W_landing_fore= Wtotal_tons + W_owtank_remaining_landing + W_iwtank_remaining_landing + W_ftank + Wpax + Wluggage;
CG_landing_fore= (Wtotal_tons * CGtotal + W_owtank_remaining_landing*CG_owtank + W_iwtank_remaining_landing*CG_iwtank + W_ftank*CG_ftank + Wpax * CGpax + Wluggage * CGluggage) / W_landing_fore;
disp(['Total CG (landing fore): ', num2str(CG_landing_fore), ' m']);

W_iwtank_left_landing =0;
W_ftank_left_landing = abs(W_iwtank - (Wfuel - (0.9700* 0.9850 *0.6225* 0.9900*0.9850* 0.9873*0.9844* 0.9900* 0.9950)*Wfuel));
W_landing_aft= Wtotal_tons + W_owtank + W_iwtank_left_landing + W_ftank_left_landing + Wpax + Wluggage;
CG_landing_aft= (Wtotal_tons * CGtotal + W_owtank*CG_owtank + W_iwtank_left_landing*CG_iwtank + W_ftank_left_landing*CG_ftank + Wpax * CGpax + Wluggage * CGluggage) / W_landing_aft;
disp(['Total CG (landing aft): ', num2str(CG_landing_aft), ' m']);

%%without passengers
W_owtank_remaining_landing=0;
W_iwtank_remaining_landing= abs(W_owtank - (Wfuel - (0.9700* 0.9850 *0.6225* 0.9900*0.9850* 0.9873*0.9844* 0.9900* 0.9950)*Wfuel)) ;
W_landing_fore_nopax= Wtotal_tons + W_owtank_remaining_landing + W_iwtank_remaining_landing + W_ftank;
CG_landing_fore_nopax= (Wtotal_tons * CGtotal + W_owtank_remaining_landing*CG_owtank + W_iwtank_remaining_landing*CG_iwtank + W_ftank*CG_ftank) / W_landing_fore_nopax;
disp(['Total CG (landing fore no pax): ', num2str(CG_landing_fore_nopax), ' m']);

W_iwtank_left_landing =0;
W_ftank_left_landing = abs(W_iwtank - (Wfuel - (0.9700* 0.9850 *0.6225* 0.9900*0.9850* 0.9873*0.9844* 0.9900* 0.9950)*Wfuel));
W_landing_aft_nopax= Wtotal_tons + W_owtank + W_iwtank_left_landing + W_ftank_left_landing;
CG_landing_aft_nopax= (Wtotal_tons * CGtotal + W_owtank*CG_owtank + W_iwtank_left_landing*CG_iwtank + W_ftank_left_landing*CG_ftank) / W_landing_aft_nopax;
disp(['Total CG (landing aft no pax): ', num2str(CG_landing_aft_nopax), ' m']);

W_fuel_nopass= Wtotal_tons + Wfuel + W_trimtank;
CG_fuel_nopass = (Wtotal_tons * CGtotal + (Wfuel + W_trimtank) * CGfuel) / W_fuel_nopass;

% figure
% plot(W_trimtank, abs(CGfuel-CGtotal_nofuel), LineWidth=2)
% xlabel('trim tank mass (Kg)')
% ylabel('difference between cg with no fuel and cg full')

LE_mac = CGw - (1/4)*mac
figure
% Set marker size
marker_size = 100;
line_width = 2;

% Plot scatter points with unique bright colors
cruise_end_cg =[((CG_cruise_end_aft_nopax-LE_mac)/mac)*100,((CG_cruise_end_fore_nopax-LE_mac)/mac)*100,((CG_cruise_end_fore-LE_mac)/mac)*100, ((CG_cruise_end_aft-LE_mac)/mac)*100, ((CG_cruise_end_aft_nopax-LE_mac)/mac)*100];
weights_cruise_end = [ W_cruise_end_aft_nopax,W_cruise_end_fore_nopax,W_cruise_end_fore,   W_cruise_end_aft,  W_cruise_end_aft_nopax];

cruise_start_cg = [((CG_cruise_start_aft_nopax-LE_mac)/mac)*100, ((CG_cruise_start_fore_nopax-LE_mac)/mac)*100, ((CG_cruise_start_fore-LE_mac)/mac)*100, ((CG_cruise_start_aft-LE_mac)/mac)*100, ((CG_cruise_start_aft_nopax-LE_mac)/mac)*100]; 
weights_cruise_start= [W_cruise_start_aft_nopax, W_cruise_start_fore_nopax,  W_cruise_start_fore, W_cruise_start_aft,W_cruise_start_aft_nopax ];

overall_envelope_cg =[((CGtotal-LE_mac)/mac)*100, ((CG_cruise_end_fore-LE_mac)/mac)*100, ((CGtotal_full-LE_mac)/mac)*100, ((CG_cruise_end_aft_nopax-LE_mac)/mac)*100, ((CGtotal-LE_mac)/mac)*100];
weight_overall_envelope = [Wtotal_tons, W_cruise_end_fore, Wtotal_full, W_cruise_end_aft_nopax, Wtotal_tons];

x_positions_landing = [
    ((CG_landing_fore - LE_mac) / mac) * 100,        % Forward CG with passengers
    ((CG_landing_aft - LE_mac) / mac) * 100,        % Aft CG with passengers
    ((CG_landing_aft_nopax - LE_mac) / mac) * 100,  % Aft CG without passengers
    ((CG_landing_fore_nopax - LE_mac) / mac) * 100, % Forward CG without passengers
    ((CG_landing_fore - LE_mac) / mac) * 100        % Close the quadrilateral
];

y_weights_landing = [
    W_landing_fore,         % Weight at forward CG with passengers
    W_landing_aft,          % Weight at aft CG with passengers
    W_landing_aft_nopax,    % Weight at aft CG without passengers
    W_landing_fore_nopax,   % Weight at forward CG without passengers
    W_landing_fore          % Close the quadrilateral
];

% Define consistent colors for the connecting lines
color_cruise_start = [0.2, 0.5, 0.8]; % Example color for cruise start line
color_cruise_end = [0.3, 0.7, 0.3];   % Example color for cruise end line
color_landing = [0.9, 0.4, 0.1];      % Example color for landing line

figure
% Empty weight
scatter(((CGtotal-LE_mac)/mac)*100, Wtotal_tons, marker_size, 'r', 'x', 'LineWidth', line_width); % Empty weight
hold on;

% Fuel, no passengers
% scatter(((CG_fuel_nopass-LE_mac)/mac)*100, W_fuel_nopass, marker_size, 'g', 'x', 'LineWidth', line_width); % Fuel, no passengers %takeoff, no pax

% Full CG (fuel and passengers)
scatter(((CGtotal_full-LE_mac)/mac)*100, Wtotal_full, marker_size, "black", 'x', 'LineWidth', line_width); % Full CG, takeoff with passengers

% Aft and front CG limits
xline(((CG_cruise_end_aft_nopax-LE_mac)/mac)*100, '--', 'Aftmost CG ','LabelHorizontalAlignment', 'left', LineWidth=1.5); % Aft CG limit
xline(((CG_cruise_end_fore-LE_mac)/mac)*100, '--', 'Foremost CG ', LineWidth=1.5); % Front CG limit

% Cruise Start Points (matching line color)
scatter(((CG_cruise_start-LE_mac)/mac)*100, W_cruise_start, marker_size, color_cruise_start, 'x', 'LineWidth', line_width); % Cruise start
% scatter(((CG_cruise_start_fore-LE_mac)/mac)*100, W_cruise_start_fore, marker_size, color_cruise_start, 'x', 'LineWidth', line_width); % Cruise start forward
% scatter(((CG_cruise_start_fore_nopax-LE_mac)/mac)*100, W_cruise_start_fore_nopax, marker_size, color_cruise_start, 'x', 'LineWidth', line_width); % Cruise start forward, no pax
% scatter(((CG_cruise_start_aft-LE_mac)/mac)*100, W_cruise_start_aft, marker_size, color_cruise_start, 'x', 'LineWidth', line_width); % Cruise start aft
% scatter(((CG_cruise_start_aft_nopax-LE_mac)/mac)*100, W_cruise_start_aft_nopax, marker_size, color_cruise_start, 'x', 'LineWidth', line_width); % Cruise start aft, no pax
% plot(cruise_start_cg, weights_cruise_start,'color',color_cruise_start); % Cruise start connecting line

% Cruise End Points (matching line color)
scatter(((CG_cruise_end-LE_mac)/mac)*100, W_cruise_end, marker_size, "blue", '^', 'LineWidth', line_width); % Cruise end
scatter(((CG_cruise_end_fore-LE_mac)/mac)*100, W_cruise_end_fore, marker_size, color_cruise_end, '^', 'LineWidth', line_width); % Cruise end forward
% scatter(((CG_cruise_end_fore_nopax-LE_mac)/mac)*100, W_cruise_end_fore_nopax, marker_size, color_cruise_end, '^', 'LineWidth', line_width); % Cruise end forward, no pax
% scatter(((CG_cruise_end_aft-LE_mac)/mac)*100, W_cruise_end_aft, marker_size, color_cruise_end, '^', 'LineWidth', line_width); % Cruise end aft
scatter(((CG_cruise_end_aft_nopax-LE_mac)/mac)*100, W_cruise_end_aft_nopax, marker_size, "black", '^', 'LineWidth', line_width); % Cruise end aft, no pax
% plot(cruise_end_cg, weights_cruise_end, 'color',color_cruise_end'); % Cruise end connecting line

% Landing Points (matching line color)
% scatter(((CG_landing_fore - LE_mac) / mac) * 100, W_landing_fore, marker_size, color_landing, 's', 'LineWidth', line_width); % Landing fore with passengers
% scatter(((CG_landing_aft - LE_mac) / mac) * 100, W_landing_aft, marker_size, color_landing, 's', 'LineWidth', line_width); % Landing aft with passengers
% scatter(((CG_landing_fore_nopax - LE_mac) / mac) * 100, W_landing_fore_nopax, marker_size, color_landing, 's', 'LineWidth', line_width); % Landing fore without passengers
% scatter(((CG_landing_aft_nopax - LE_mac) / mac) * 100, W_landing_aft_nopax, marker_size, color_landing, 's', 'LineWidth', line_width); % Landing aft without passengers
scatter(((CG_landing_design - LE_mac) / mac) * 100, W_landing, marker_size, color_landing, 's', 'LineWidth', line_width);
% plot(x_positions_landing, y_weights_landing, 'color',color_landing'); % Landing connecting line

fill(overall_envelope_cg, weight_overall_envelope, 'cyan', 'FaceAlpha','0.1')

% Add weight lines
yline(Wtotal_full + W_trimtank, '--k', 'MTOW', 'LabelHorizontalAlignment','center', LineWidth=1.5); % Maximum takeoff weight
yline(Wtotal_tons, '--k', 'Empty Weight','LabelHorizontalAlignment','center', LineWidth=1.5); % Empty weight

xline(((CGtotal_full-LE_mac)/mac)*100, '--m', 'Operational condtions',LineWidth=1.5)
xline(((33.9- LE_mac)/mac)*100, '-r', 'Forward cg limit', LineWidth=2)
xline(((35.32- LE_mac)/mac)*100, '-r', 'Aft cg limit', LineWidth=2)

% Add legend
legend({
    'Empty weight', ...
    'Full CG (Fuel and Pax)', ...
    'Aftmost CG', ...
    'Foremost CG', ...
    'Cruise Start', ...
    'Cruise End', ...
    'Cruise End Forward', ...
    'Cruise End Aft (No Pax)', ...
    'Landing', 'CG envelope'...
    'MTOW (Max Takeoff Weight)', ...
    'Empty Weight', ...
    'Operational conditions', 'Forward cg limit', 'Aft cg limit'}, ...
    'Location', 'bestoutside');
grid on
ylabel('Mass (tonnes)')
xlabel('CG position (% MAC)')

% figure
% % scatter(((CGtotal-LE_mac)/mac)*100, Wtotal_tons,'cyan') %empty weight
% hold on
% scatter(CG_cruise_end, W_cruise_end, 'blue')
% scatter(CGtotal_nofuel, Wtotal_nofuel, 'red')
% scatter(CGtotal_full, Wtotal_full, 'magenta')
% legend("Cruise Start", "Cruise End", "0% Fuel", "100% Fuel")


% get Iyy using sum(Iyy) = sum(m * (x - x_cg)^2)
Iyy = Ww * (CGw - CGtotal)^2 + Wht * (CGht - CGtotal)^2 + Wvt * (CGvt - CGtotal)^2 + ...
    Wfus * (CGfus - CGtotal)^2 + Wmlg * (CGmlg - CGtotal)^2 + Wnlg * (CGnlg - CGtotal)^2 + ...
    Winl * (CGinl - CGtotal)^2 + Wec * (CGec - CGtotal)^2 + Wes * (CGes - CGtotal)^2 + ...
    Wfs * (CGfs - CGtotal)^2 + Wfc * (CGfc - CGtotal)^2 + WAPUinst * (CGapu - CGtotal)^2 + ...
    Winstr * (CGinstr - CGtotal)^2 + Whydr * (CGhydr - CGtotal)^2 + Wel * (CGel - CGtotal)^2 + ...
    Wav * (CGav - CGtotal)^2 + Wfurn * (CGfurn - CGtotal)^2 + Wac * (CGac - CGtotal)^2 + ...
    Wai * (CGai - CGtotal)^2 + Whg_civilian * (CGwhg - CGtotal)^2 + Whg_military * (CGmil - CGtotal)^2 + ...
    W_engines * (CGeng - CGtotal)^2 + Wfuel * (CGfuel - CGtotal_full)^2 + ...
    Wpax * (CGpax - CGtotal_full)^2 + Wluggage * (CGluggage - CGtotal_full)^2 + ...
    (1/(12 * 3.208^2)) * Wfus * (1.5 * D^2 + L^2);

% convert from lb*m^2 to kg*m^2
Iyy = Iyy * 0.453592;

disp(['Iyy: ', num2str(Iyy), ' kg*m^2']);

%get the data for the table

%%takeoff
%with passesngers

takeoff_fore = ((CGtotal_full-LE_mac)/mac)*100;
cruise_start_fore_pax = ((CG_cruise_start_fore-LE_mac)/mac)*100;
cruise_start_aft_pax= ((CG_cruise_start_aft-LE_mac)/mac)*100;
cruise_end_fore_pax = ((CG_cruise_end_fore-LE_mac)/mac)*100;
cruise_end_aft_pax = ((CG_cruise_end_aft-LE_mac)/mac)*100;
landing_fore = ((CG_landing_fore-LE_mac)/mac)*100;
landing_aft= ((CG_landing_aft-LE_mac)/mac)*100;

%without passengers

takeoff_fore_nopax = ((CG_fuel_nopass-LE_mac)/mac)*100;
cruise_start_fore_nopax = ((CG_cruise_start_fore_nopax-LE_mac)/mac)*100;
cruise_start_aft_nopax= ((CG_cruise_start_aft_nopax-LE_mac)/mac)*100;
cruise_end_fore_nopax = ((CG_cruise_end_fore_nopax-LE_mac)/mac)*100;
cruise_end_aft_nopax = ((CG_cruise_end_aft_nopax-LE_mac)/mac)*100;
landing_fore_nopax = ((CG_landing_fore_nopax-LE_mac)/mac)*100;
landing_aft_nopax= ((CG_landing_aft_nopax-LE_mac)/mac)*100;

%% 
% % Constants and aircraft configuration
% empty_weight_CG = CGtotal; % CG position of the empty aircraft (meters)
% empty_weight = Wtotal_tons; % Empty operating weight of the aircraft (tons)
% 
% % Fuel tank properties
% tank_data = {
%     % Tank Name      Position_CG (m)   Max_Fuel_Weight (tons)
%     'Owtank     ',   CG_owtank,        W_owtank; 
%     'Innerwing',   CG_iwtank,            W_iwtank; 
%     'Center Tank',   CG_ftank,             W_ftank ; 
% };
% 
% % Extract data for calculations
% tank_positions = cell2mat(tank_data(:, 2)); % Tank CG positions
% tank_capacities = cell2mat(tank_data(:, 3)); % Tank max fuel weights
% total_fuel_capacity = sum(tank_capacities)*0.97*0.985; % Total fuel capacity
% 
% % Initialize outputs for foremost and aftmost CG calculations
% foremost_CG = []; % Forward CG limit
% aftmost_CG = []; % Aft CG limit
% 
% % Loop through symmetric fuel burn scenarios
% fuel_burn_percentages = linspace(0, 1, 50); % From full fuel (0%) to empty (100%)
% for i = 1:length(fuel_burn_percentages)
%     fuel_fraction = 1 - fuel_burn_percentages(i); % Remaining fuel fraction
% 
%     % Symmetric fuel burn: Reduce fuel evenly across all tanks
%     fuel_weights = fuel_fraction * tank_capacities;
% 
%     % Calculate total weight and CG
%     total_weight = empty_weight + sum(fuel_weights); % Aircraft weight (tons)
%     weighted_CG = (empty_weight_CG * empty_weight + sum(fuel_weights .* tank_positions)) / total_weight; % CG position (meters)
% 
%     % Append CG for current fuel distribution
%     foremost_CG = [foremost_CG, weighted_CG]; % Foremost CG progression
%     aftmost_CG = [aftmost_CG, weighted_CG]; % Aftmost CG progression
% end
% 
% % Find extreme CG positions
% foremost_CG_position = min(foremost_CG); % Foremost CG (smallest CG value)
% aftmost_CG_position = max(aftmost_CG); % Aftmost CG (largest CG value)
% 
% % Convert CG positions to % MAC
% foremost_CG_percent_MAC = ((foremost_CG_position - LE_mac) / mac) * 100;
% aftmost_CG_percent_MAC = ((aftmost_CG_position - LE_mac) /mac) * 100;
% 
% % Display results
% disp(['Foremost CG: ', num2str(foremost_CG_percent_MAC), '% MAC']);
% disp(['Aftmost CG: ', num2str(aftmost_CG_percent_MAC), '% MAC']);
% 
% % Plot CG variation over fuel burn
% figure;
% plot(fuel_burn_percentages * 100, ((foremost_CG - LE_mac) / mac) * 100, 'b', 'LineWidth', 2);
% hold on;
% plot(fuel_burn_percentages * 100, ((aftmost_CG - LE_mac) / mac) * 100, 'r', 'LineWidth', 2);
% xlabel('Fuel Burn Percentage (%)');
% ylabel('CG Position (% MAC)');
% title('CG Position vs. Fuel Burn');
% legend('Foremost CG', 'Aftmost CG');
% grid on;
% hold off;
% 
% 
