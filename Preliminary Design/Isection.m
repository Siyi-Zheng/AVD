function [A , I] = Isection(t , b , h)
    
    %calculate the required flange thickness​
    A = 2 * t * b + (h - 2 *t) * t;

    %check moment of inertia
    I = (A * 3.17^2) / 2;

end