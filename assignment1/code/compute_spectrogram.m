function [stft_mat, freq_axis, time_axis] = compute_spectrogram(x, Fs, win_len, win_shift)
    x_len = size(x,1);
    num_seg = ceil((x_len - win_len)/win_shift);
    fft_len = win_len;

    stft_mat = zeros(fft_len, num_seg);
%     spec_double_sided = zeros(fft_len, num_seg);
%     spec_single_sided = zeros(fft_len/2 + 1, num_seg);

    hamm_win_samp = hamming(win_len);
    
    % iterate over all segments and compute dft for each segment
    for r = 0:num_seg-1
        win_start = r * win_shift + 1;
        win_end = r * win_shift + win_len;
        
        if win_end > x_len
            x_seg = [x(win_start:x_len); zeros(win_end-x_len, 1)] .* hamm_win_samp;
        else
            x_seg = x(win_start : win_end) .* hamm_win_samp;
        end
        
        assert(size(x_seg, 1) == win_len);
        x_seg_dft = fft(x_seg);
%         x_seg_dft_len = size(x_seg_dft, 1);

        % compute the single-sided spectrum
        stft_mat(:, r+1) = x_seg_dft;
%         spec_double_sided(:, r+1) = 20 * log10(abs(x_seg_dft));
%         spec_single_sided(:, r+1) = spec_double_sided(1:x_seg_dft_len/2 + 1, r+1);
    end

    time_axis = linspace(0, num_seg-1) * win_shift / Fs;
    freq_axis = (Fs * (0:fft_len/2)) / (fft_len * 1000);
end