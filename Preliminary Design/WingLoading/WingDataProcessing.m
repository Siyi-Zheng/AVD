clear 
clc

%considering the effect of what altitude the analysis is conducted at for
%wing loading

heightsm = 0:100:40000;
H = convlength(heightsm,'ft','m');

Sref = 488.545;

format long g;

[T,a,P,rhos,nu] = atmosisa(H);

mu = nu .* rhos;

Vd = convvel(333 , 'kts' , 'm/s');
Va = convvel(254 , 'kts' , 'm/s');

Mach_noa = Va ./ a;
Mach_nod = Vd ./a;

Rea = (12.41 .* rhos .* Va) ./ mu; 
Red = (12.41 .* rhos .* Vd) ./ mu;

%lift for 2.5g turn

L = 354000 * 2.5 * 9.81;

Cl_requireda = L ./ (0.5 * Sref * Va^2 .* rhos);
Cl_requiredd = L ./ (0.5 * Sref * Vd^2 .* rhos);




%%

%wing load data processing

%define air parameters
v = 130.7;
rho = 1.225;

%load AOA
alpha = readmatrix("LoadCase1.csv" , 'Range', 'A12:A68')';

%load wing data from CSV file
Spans = readmatrix("MainWing_a=-4.00_v=130.70ms.csv" , 'Range' , 'A22:A265')';
Chords = readmatrix("MainWing_a=-4.00_v=130.70ms.csv" , 'Range' , 'B22:B265')';
Cl(1,:) = readmatrix("MainWing_a=-4.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(2,:) = readmatrix("MainWing_a=-3.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(3,:) = readmatrix("MainWing_a=-3.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(4,:) = readmatrix("MainWing_a=-2.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(5,:) = readmatrix("MainWing_a=-2.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(6,:) = readmatrix("MainWing_a=-1.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(7,:) = readmatrix("MainWing_a=-1.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(8,:) = readmatrix("MainWing_a=-0.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(9,:) = readmatrix("MainWing_a=0.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(10,:) = readmatrix("MainWing_a=0.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(11,:) = readmatrix("MainWing_a=1.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(12,:) = readmatrix("MainWing_a=1.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(13,:) = readmatrix("MainWing_a=2.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(14,:) = readmatrix("MainWing_a=2.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(15,:) = readmatrix("MainWing_a=3.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(16,:) = readmatrix("MainWing_a=3.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(17,:) = readmatrix("MainWing_a=4.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(18,:) = readmatrix("MainWing_a=4.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(19,:) = readmatrix("MainWing_a=5.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(20,:) = readmatrix("MainWing_a=5.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(21,:) = readmatrix("MainWing_a=6.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(22,:) = readmatrix("MainWing_a=6.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(23,:) = readmatrix("MainWing_a=7.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(24,:) = readmatrix("MainWing_a=7.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(25,:) = readmatrix("MainWing_a=8.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(26,:) = readmatrix("MainWing_a=8.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(27,:) = readmatrix("MainWing_a=9.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(28,:) = readmatrix("MainWing_a=9.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(29,:) = readmatrix("MainWing_a=10.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(30,:) = readmatrix("MainWing_a=10.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(31,:) = readmatrix("MainWing_a=11.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(32,:) = readmatrix("MainWing_a=11.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(33,:) = readmatrix("MainWing_a=12.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(34,:) = readmatrix("MainWing_a=12.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(35,:) = readmatrix("MainWing_a=13.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(36,:) = readmatrix("MainWing_a=13.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(37,:) = readmatrix("MainWing_a=14.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(38,:) = readmatrix("MainWing_a=14.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(39,:) = readmatrix("MainWing_a=15.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(40,:) = readmatrix("MainWing_a=15.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(41,:) = readmatrix("MainWing_a=16.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(42,:) = readmatrix("MainWing_a=16.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(43,:) = readmatrix("MainWing_a=17.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(44,:) = readmatrix("MainWing_a=17.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(45,:) = readmatrix("MainWing_a=18.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(46,:) = readmatrix("MainWing_a=18.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(47,:) = readmatrix("MainWing_a=19.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(48,:) = readmatrix("MainWing_a=19.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(49,:) = readmatrix("MainWing_a=20.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(50,:) = readmatrix("MainWing_a=20.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(51,:) = readmatrix("MainWing_a=21.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(52,:) = readmatrix("MainWing_a=21.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(53,:) = readmatrix("MainWing_a=22.00_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(54,:) = readmatrix("MainWing_a=22.50_v=130.70ms.csv" , 'Range' , 'D22:D265');
Cl(55,:) = readmatrix("MainWing_a=23.00_v=130.70ms.csv" , 'Range' , 'D22:D265');




%removing fuselage lift
I = (-2.55 < Spans) & (Spans < 2.55);
I = find(I == 1);
Cl(: , I) = 0;


%removing lift in central section and considering only one half of the wing
I = find(Spans>2.55, 1 );
Cl_half = Cl( : , 123 : end);

% = 123;

%converting to sectional lift
for i = 1:55

    Sect_L(i,:) = Cl(i,:) * 0.5 * rho * v^2 .* Chords;
    Sect_Lhalf(i,:) = Cl_half(i,:) * 0.5 * rho * v^2 .* Chords(123:end);

end


%plotting
figure
hold on

for i = 1:size(Sect_L , 1)

    plot(Spans , Sect_L , LineWidth=2)

end

hold off

%need to find which lift dist represents a lift of 2.5g

%finding the lifts of each aoa

Lifts = zeros(1, 55); % Where `n` is the number of iterations

for i = 1:size(Sect_Lhalf , 1)

    Lifts(i) = 2 * trapz(Spans(123:end)' , Sect_Lhalf(i , :));

end

Cls = (Lifts ./ (0.5 * rho * v^2 * Sref))';

%best AoA to abtain a loading of 2.5g
[~ , ind] = min(abs(Lifts - L));

aoa_opt = alpha(ind);

disp(['The aoa to get the correct load factor for Va is ' , num2str(aoa_opt) , ' degrees.'])

%Converting to correct disctretisation
xq = [2.55 : 0.01 : 32.5 , 32.5];

vq = interp1(Spans(123:end) , Sect_Lhalf(ind,:) , xq , "pchip");

%vq(1:255) = 0;

figure
hold on

plot(Spans(123:end) , Sect_Lhalf(ind , :) , 'LineStyle','none' , 'Marker','*' ,  'DisplayName', 'Data')
plot(xq , vq , 'LineStyle','--' , 'Marker' , 'none' ,  'DisplayName', 'Fit')


hold off

vq = vq .* 1.0228;

%check it's correct
Lift_fit = 2 * trapz(xq , vq);

save("WingLoadCase1.mat" , "vq" );

%scaled wing load for 1g
vq = vq*(1/2.5);

%check it's correct
Lift_fit = 2 * trapz(xq , vq);

save("WingLoadCase2.mat" , "vq");


