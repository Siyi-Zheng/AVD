function [M_new , N_new , S_new] = anglechange(M , N , S , angle)

    index = 3601 - (angle / 0.1 + 1);

    M_new = [M((index):end) , M(1:(index-1))];
    N_new = [N((index):end) , N(1:(index-1))];
    S_new = [S((index):end) , S(1:(index-1))];

end