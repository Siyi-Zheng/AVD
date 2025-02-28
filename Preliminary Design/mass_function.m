function mass = mass_function(no_stringer, A_stringer, Lfs, t_skin)
    % Constants
    l_plane = 77.8;         
    d_fuslg = 6.34;  
    circum_fuslg = pi*d_fuslg;

    % Compute Af_min from fuslg_failure
    [~, ~, ~, Af_min] = fuslg_failure(no_stringer, A_stringer, Lfs, t_skin);

    % Compute mass
    mass = 2765*l_plane * t_skin * circum_fuslg + ...
           2790*A_stringer * no_stringer * l_plane + ...
           2765*Af_min * circum_fuslg * (l_plane / Lfs);
end
