skin = t_skin*pi*6.3*77.8*2765

stringer = A_stringer*77.8*no_stringer*2765

[~, ~, ~, Af_min] = fuslg_failure(no_stringer, A_stringer, Lfs, t_skin);


light_frame = Af_min*pi*6.3*(77.8/Lfs)*2765
