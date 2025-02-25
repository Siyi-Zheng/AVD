function [c, ceq] = fuslg_constraints(x)
    % Extract optimization variables
    no_stringer = round(x(1));
    A_stringer = x(2);
    Lfs = x(3);

    % Constants
    sig_yield = 431e6;

    % Get values from fuslg_failure
    [sig_skin, sig_euler, sig_max, ~] = fuslg_failure(no_stringer, A_stringer, Lfs);

    % Define inequality constraints (should be ≤ 0)
    c = [
        sig_max - sig_yield;  % sig_max < sig_yield
        sig_max - sig_skin;   % sig_max < sig_skin
        sig_max - sig_euler   % sig_max < sig_euler
    ];

    % No equality constraints
    ceq = [];
end
