function [mlf , If , nf , Af , t] = lightframes(Lfs,h,b,sectionshape,E,rho)
%function to calculate total frame mass for a given frame spacing and
%dimensions.
%first compute frame inertia, Ilf such that panel and general instability
%occur under the same loading conditions, provided

D = 6.34; %frame diameter

Cf = 1/16000; % empirical based off work performed by Lockheed Martin and Shanley, J. Aero Science

%syms Lfs %frame spacing

M_ult = 1.6e7;

%set a fixed IF_xx
If = (Cf * M_ult * D^2) / (E * Lfs);


%C design or rectangular design (delete as needed)
%syms t;
%syms b;
%syms h;

if sectionshape == true
    b = ((If * 12) / h^3);
    Af = b * h;
else
    t = If / ((h^3/12) + (0.5 * b * h^2));
    Af = (2 * b + h) * t;
end

nf = 51.6/Lfs;

%calculating total frame mass
mlf = rho * (Af * pi * D) * nf;

end
