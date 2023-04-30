function [x_reconst] = compute_reconst_sig(stft_mat, x_len, win_len, win_shift)

    [fft_len, num_seg] = size(stft_mat);
    x_hat = zeros(x_len, 1);
%     stft_mag = 10.^(spec / 20.0);
    
    % iterate over each DFT segment to reconstruct the original sequence
    for r = 0:num_seg-1
        x_seg_dft = stft_mat(:, r+1);

%         x_seg_dft_mag_double = zeros(2 * (fft_len-1), 1);
%         x_seg_dft_mag_double(1:fft_len-1) = x_seg_dft_mag(2:end);
%         x_seg_dft_mag_double(fft_len:end) = flip(conj(x_seg_dft_mag(2:end)));
%         x_seg_dft_mag = spec(:, r+1);
%         x_seg = ifft(x_seg_dft_mag_double);

        x_seg = ifft(x_seg_dft);

        win_start = r * win_shift + 1;
        win_end = r * win_shift + win_len;

        if win_end > x_len
            x_hat(win_start:x_len) = x_hat(win_start:x_len) + x_seg(1:x_len - win_start + 1);
        else
            x_hat(win_start:win_end) = x_hat(win_start:win_end) + x_seg;
        end
    end

    % construct a repeated version of weight function (since it is periodic)
    weight_func = get_weight_func(win_len, win_shift, "hamming");
    weight_func_repeat = zeros(x_len, 1);
    for i = 1:win_len:x_len
        if i + win_len - 1 > x_len
            weight_func_repeat(i:x_len) = weight_func(1:x_len - i + 1);
        else
            weight_func_repeat(i:i + win_len - 1) = weight_func;
        end
    end

    assert(size(x_hat, 1) == size(weight_func_repeat, 1));
    x_reconst = x_hat ./ weight_func_repeat;
end