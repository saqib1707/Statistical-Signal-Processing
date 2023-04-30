function [corr_mat] = compute_corrmat(ryy, P)
    N = (size(ryy, 1) + 1) / 2;

    corr_mat = zeros(P, P);
    for i = 1:P
        corr_mat(i, :) = ryy(N - i + 1 : N - i + P);
    end

    first_row = ryy(N:N+P-1);
    corr_mat2 = toeplitz(first_row);

    assert(size(corr_mat,1) == size(corr_mat2,1));
    assert(size(corr_mat,2) == size(corr_mat2,2));
    assert( sum(abs(corr_mat2 - corr_mat), 'all') < 1e-10 );

end