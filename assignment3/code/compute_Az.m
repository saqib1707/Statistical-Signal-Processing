function [Az] = compute_Az(ap_par_vec, P, N)
    assert(size(ap_par_vec, 1) == P+1);
    Az = zeros(N, 1);
    l_vec = transpose(0:P);
    
    for n = 0:N-1
        omega = 2 * pi * n / N;
        tmp1 = conj(ap_par_vec);
        tmp2 = exp(-1j * omega .* l_vec);
%         assert(size(tmp1) == size(tmp2));
        Az(n+1, 1) = sum(tmp1 .* tmp2);
    end
end