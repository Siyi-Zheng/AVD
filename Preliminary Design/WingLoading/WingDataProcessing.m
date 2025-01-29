clear 
clc

%considering the effect of what altitude the analysis is conducted at for
%wing loading

heightsm = 0:100:40000;
H = convlength(heightsm,'ft','m');

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

Cl_requireda = L ./ (0.5 * 488.545 * Va^2 .* rhos);
Cl_requiredd = L ./ (0.5 * 488.545 * Vd^2 .* rhos);




%%

%wing load data processing

%define air parameters
v = 130.7;
rho = 1.225;

%load AOA
alpha = readmatrix("LoadCase1.csv" , 'Range', 'A8:A46');

%load wing data from CSV file
Spans = readmatrix("MainWing_a=-4.00_v=130.70ms.csv" , 'Range' , 'A22:A265');
Chords = readmatrix("MainWing_a=-4.00_v=130.70ms.csv" , 'Range' , 'B22:B265');
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

%removing fuselage lift
I = (-2.55 < Spans) & (Spans < 2.55);
I = find(I == 1);
Cl(: , I) = 0;


%removing lift in central section and considering only one half of the wing
I = find(Spans>2.55, 1 );
Cl_half = Cl( : , I : end);

%converting to sectional lift
for i = 1:39

    Sect_L(i,:) = Cl(i,:) * 0.5 * rho * v^2 .* Chords';
    Sect_Lhalf(i,:) = Cl_half(i,:) * 0.5 * rho * v^2 .* Chords(I:end)';

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

Lifts = zeros(1, 39); % Where `n` is the number of iterations

for i = 1:size(Sect_Lhalf , 1)

    Lifts(i) = 2 * trapz(Sect_Lhalf(i , :));

end



%best AoA to abtain a loading of 2.5g
[~ , ind] = min(abs(Lifts - L));

aoa_opt = alpha(ind);

disp(['The aoa to get the correct load factor for Va is ' , num2str(aoa_opt) , ' degrees.'])

%Converting to correct disctretisation
xq = 2.55 : 0.01 : 32.5;

vq = interp1(Spans(I:end) , Sect_Lhalf(ind,:) , xq , "pchip");

%vq(1:255) = 0;

figure
hold on

plot(Spans(I:end) , Sect_Lhalf(ind , :) , 'LineStyle','none' , 'Marker','*' ,  'DisplayName', 'Data')
plot(xq , vq , 'LineStyle','--' , 'Marker' , 'none' ,  'DisplayName', 'Fit')

xlim([0 40])
ylim([-20000 160000])

hold off

%check it's correct
Lift_fit = 2 * trapz(xq , vq);

data = Sect_Lhalf(ind , :);

save("WingLoad.mat" , "data" )



