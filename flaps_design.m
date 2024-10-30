clear
clc
c_prime = 1.0448; % normalised to chord length
te_kink = 0.3; % as a proportion of the span
sweep = 30;
S_ratio = 0.75; % proportion of wing area affected by flaps

delta_cl_inner = 0.9 * 1.6 * c_prime;
delta_cl_outer = 0.9 * S_ratio * cosd(sweep) * 1.6 * c_prime;
delta_cl = delta_cl_inner * te_kink + delta_cl_outer * (1 - te_kink);