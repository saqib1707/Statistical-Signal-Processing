function [autocorr_mat] = compute_autocorr_toeplitz_mat(ryy, P)
    N = (size(ryy, 1) + 1) / 2;

    autocorr_mat = zeros(P+1, P+1);
    for i = 1:P+1
        autocorr_mat(i, :) = ryy(N - i + 1 : N - i + 1 + P);
    end

    first_row = ryy(N:N+P);
    autocorr_mat2 = toeplitz(first_row);

    assert(size(autocorr_mat,1) == size(autocorr_mat2,1));
    assert(size(autocorr_mat,2) == size(autocorr_mat2,2));
    assert( sum(abs(autocorr_mat2 - autocorr_mat), 'all') == 0 );

end