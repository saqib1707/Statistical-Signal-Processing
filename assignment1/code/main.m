clc; clear;

run_q12 = false;
run_q3 = false;
run_q4 = false;
run_q5 = true;

if run_q12 == true || run_q3 == true || run_q4 == true
    [x_p1, Fs] = audioread("part1.wav");
    
    x_p1_len = size(x_p1, 1);
    x_p1_dur = x_p1_len / Fs;
    
    % play the signal
    % sound(x_p1, Fs);
end

% -------------------------------------------------------------------------
if run_q12 == true
    tosave = true;
    for win_len = [256, 1024, 4096]
        win_shift = win_len / 2;
    
        % compute the spectrogram using in-built spectrogram function
        num_fft_samp = max(256, 2^nextpow2(win_len));
        [stft_inbuilt, freq_axis, time_axis] = spectrogram(x_p1, hamming(win_len), win_shift, num_fft_samp, Fs, 'yaxis');
        spec_single_sided_inbuilt = 20 * log10(abs(stft_inbuilt));
        plot_spectrogram(spec_single_sided_inbuilt, freq_axis/1000, time_axis, win_len, win_shift, tosave, "1_2", "spec_inbuilt_M"+win_len+"_R"+win_shift+".png");
    
        % compute the spectrogram using implemented function
        [stft_own, freq_axis, time_axis] = compute_spectrogram(x_p1, Fs, win_len, win_shift);
        [fft_len, num_seg] = size(stft_own);

        spec_double_sided_own = 20 * log10(abs(stft_own));
        spec_single_sided_own = spec_double_sided_own(1:fft_len/2 + 1, :);

        plot_spectrogram(spec_single_sided_own, freq_axis, time_axis, win_len, win_shift, tosave, "1_2", "spec_own_M"+win_len+"_R"+win_shift+".png")
    
        mse = mean((spec_single_sided_own(:) - spec_single_sided_inbuilt(:)).^2);
        disp("MSE between in-built and own spectrogram: "+mse);
    end
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
if run_q3 == true
    % plot the reconstruction weights for hamming window
    tosave = true;
    win_len = 256;
    win_shift = 128;
    weight_func = get_weight_func(win_len, win_shift, "hamming");
    plot_weight_func(weight_func, win_len, win_shift, tosave, "3", "weight_func_hamming_M"+win_len+"_R"+win_shift+".png");
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
if run_q4 == true
    % reconstruct signal from spectrogram and check whether exact recovery is achieved
    tosave = true;
    win_len = 256;
    win_shift = 128;

    num_fft_samp = max(256, 2^nextpow2(win_len));
%     [stft_mat, freq_axis, time_axis] = spectrogram(x_p1, hamming(win_len), win_shift, num_fft_samp, Fs, 'yaxis');
    [stft_mat, freq_axis, time_axis] = compute_spectrogram(x_p1, Fs, win_len, win_shift);
%     spec_single_sided = 20 * log10(abs(stft_mat));
    
    x_p1_reconst = compute_reconst_sig(stft_mat, x_p1_len, win_len, win_shift);

    mse_reconst = mean((x_p1(:) - x_p1_reconst(:)).^2);
    disp("Mean Square Error: "+mse_reconst);

    plot_reconst_sig(x_p1, x_p1_reconst, Fs, win_len, win_shift, tosave, "4");
    sound(x_p1_reconst, Fs);
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
if run_q5 == true
    [x_p5, Fs] = audioread("part5.wav");

    x_p5_len = size(x_p5, 1);
    x_p5_dur = x_p5_len / Fs;

    % sound(x_p5, Fs);

    % plot corrupted signal in time-domain
    fig = figure;
    time_axis = (1:size(x_p5, 1)) / Fs;
    plot(time_axis, x_p5);
    ylim([-1.0 1.0]);
    grid on;
    title("Corrupted Signal", "FontSize", 18);
    xlabel("Time [s]", "FontSize", 16);
    ylabel("Amplitude y[n]", "FontSize", 16);
    saveas(fig, "../plots/prob5/corrupt_sig.png");
    close;

    % compute spectrogram of corrupted signal
    win_len = 256;
    win_shift = 128;
    
    [stft_corrupt, freq_axis, time_axis] = compute_spectrogram(x_p5, Fs, win_len, win_shift);
    [fft_len, num_seg] = size(stft_corrupt);   % (256 x 1630)

    stft_mag_sq_corrupt = abs(stft_corrupt).^2;
    spec_double_sided_corrupt = 20 * log10(abs(stft_corrupt));
    spec_single_sided_corrupt = spec_double_sided_corrupt(1:fft_len/2 + 1, :);

    % plot spectrogram of corrupted signal
    plot_spectrogram(spec_single_sided_corrupt, freq_axis, time_axis, win_len, win_shift, true, "5", "spec_corrupt_M"+win_len+"_R"+win_shift+".png");

    % plot STFT magnitude squared of corrupt signal
    fig = figure;
    imagesc(time_axis, freq_axis, stft_mag_sq_corrupt(1:fft_len/2 + 1, :));
    set(gca, 'ydir', 'normal');   % flip the Y Axis so lower frequencies are at the bottom
    title("STFT Magnitude Corrupt (M = "+win_len+", R = "+win_shift+")", 'FontSize', 18);
    xlabel("Time [s]", 'FontSize', 16);
    ylabel("Frequency [KHz]", 'FontSize', 16);
    % colormap(jet(256));
    c = colorbar;
    c.Label.String = 'Intensity';
    c.Label.FontSize = 16;
    saveas(fig, "../plots/prob"+5+"/"+"stft_magsq_corrupt_M"+win_len+"_R"+win_shift+".png");
    close;

    % create a copy of corrupt STFT for processing
    stft_clean = stft_corrupt;

    % remove noise from corrupted signal
    clean_val = 0;
    stft_clean(4:6, 280:656) = clean_val;
    stft_clean(252:254, 280:656) = clean_val;
    stft_clean(1:4, 525:1240) = clean_val;
    stft_clean(254:256, 525:1240) = clean_val;

%     stft_clean(1:4, 175:240) = clean_val;
%     stft_clean(253:256, 175:240) = clean_val;
%     stft_clean(1:2, 275:350) = clean_val;
%     stft_clean(255:256, 275:350) = clean_val;
%     stft_clean(1:3, 110:130) = clean_val;
%     stft_clean(254:256, 110:130) = clean_val;
%     stft_clean(1:2, 1240:1290) = clean_val;
%     stft_clean(255:256, 1240:1290) = clean_val;
%     stft_clean(1:3, 1400:1470) = clean_val;
%     stft_clean(254:256, 1400:1470) = clean_val;

%     noise_loc = find(stft_clean > 10);
%     stft_clean(abs(stft_clean).^2 > 50) = 0;

    stft_mag_sq_clean = abs(stft_clean).^2;
    spec_double_sided_clean = 20 * log10(abs(stft_clean));
    spec_single_sided_clean = spec_double_sided_clean(1:fft_len/2 + 1, :);

    % plot spectrogram of clean signal
    plot_spectrogram(spec_single_sided_clean, freq_axis, time_axis, win_len, win_shift, true, "5", "spec_clean_M"+win_len+"_R"+win_shift+".png");

    % plot STFT magnitude squared of clean signal
    fig = figure;
    imagesc(time_axis, freq_axis, stft_mag_sq_clean(1:fft_len/2 + 1, :));
    set(gca, 'ydir', 'normal');   % flip the Y Axis so lower frequencies are at the bottom
    title("STFT Magnitude Clean (M = "+win_len+", R = "+win_shift+")", 'FontSize', 18);
    xlabel("Time [s]", 'FontSize', 16);
    ylabel("Frequency [KHz]", 'FontSize', 16);
    % colormap(jet(256));
    c = colorbar;
    c.Label.String = 'Intensity';
    c.Label.FontSize = 16;
    saveas(fig, "../plots/prob"+5+"/"+"stft_magsq_clean_M"+win_len+"_R"+win_shift+".png");
    close;

    fig = figure;
    freq_axis_double = (Fs * (0:fft_len)) / (fft_len * 1000);
    imagesc(time_axis, freq_axis_double, stft_mag_sq_clean);
    set(gca, 'ydir', 'normal');   % flip the Y Axis so lower frequencies are at the bottom
    title("STFT Magnitude Clean (M = "+win_len+", R = "+win_shift+")", 'FontSize', 18);
    xlabel("Time [s]", 'FontSize', 16);
    ylabel("Frequency [KHz]", 'FontSize', 16);
    % colormap(jet(256));
    c = colorbar;
    c.Label.String = 'Intensity';
    c.Label.FontSize = 16;
    saveas(fig, "../plots/prob"+5+"/"+"stft_double_sided_magsq_clean_M"+win_len+"_R"+win_shift+".png");
    close;

    % reconstruct clean version of corrupted signal
    x_p5_reconst = compute_reconst_sig(stft_clean, x_p5_len, win_len, win_shift);
%     x_p5_reconst = real(x_p5_reconst);
    sound(x_p5_reconst, Fs);
    plot_reconst_sig(x_p5, x_p5_reconst, Fs, win_len, win_shift, true, 5);
end
% -------------------------------------------------------------------------

