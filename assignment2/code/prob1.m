clc; clear;

% compute the impulse function h[n] from the transfer function H(z)
b = [1, -0.9, 0.81];
a = [1, -2.76, 3.809, -2.654, 0.924];

[hn, tn] = impz(b, a);
Nh = size(hn, 1);

fig = figure;
plot(tn, hn, LineWidth=2);
xlabel("index [n]", FontSize=16);
ylabel("h[n]", FontSize=16);
title("Impulse Response h[n]", FontSize=16);
grid on;
saveas(fig, "../plots/prob1a/impulse_response.png");
close;

% compute true autocorrelation of y[n]
savepath = "../plots/prob1a/true_autocorr.png";
ryy_true = compute_true_autocorr(hn, true, savepath);

% generate 1024 samples of white gaussian noise random process x[n] ~
% N(0,1)
Ny = 4097;
[xn, yn] = generate_random_process(hn, Ny, true);

% estimate autocorrelation sequence of y[n] for various sample lengths
num_samp_lst = [64, 128, 256, 512, 1024];
for num_samp = num_samp_lst
    % sample y[n] from the center
    y_samp = yn((Ny+1)/2 - num_samp/2 : (Ny+1)/2 + num_samp/2 - 1, 1);
    assert(size(y_samp, 1) == num_samp);

    savepath = "../plots/prob1b/est_autocorr_"+num_samp+".png";
    ryy_est = estimate_autocorr(y_samp, true, savepath);
end

% compute true power spectral density of y[n], Ryy(e^jw)
savepath = "../plots/prob1c/true_psd.png";
Ryy_true = compute_true_periodogram(ryy_true, true, savepath);


% compute periodogram estimates
num_samp_lst = [64, 128, 256, 512, 1024];
for num_samp = num_samp_lst
%     y_samp = yn(1:num_samp, 1);
    y_samp = yn((Ny+1)/2 - num_samp/2 : (Ny+1)/2 + num_samp/2 - 1, 1);
    assert(size(y_samp, 1) == num_samp);

    savepath = "../plots/prob1d/est_periodogram_"+num_samp+".png";
    Ryy_est = estimate_periodogram(y_samp, true, savepath);

%     ryy_est = estimate_autocorr(y_samp, false, "");
%     Ryy_est_2 = abs(fft(ryy_est));
%     close;
%     disp("This should be same:"+max(abs(Ryy_est_2 - Ryy_est)));
end


% generate several realization of the periodogram estimate
num_realizations = 1000;
Ny = 4097;
num_samp = 1024;

% estimate sample mean of the periodogram estimate
savepath = "../plots/prob1e/periodogram_mean_comparison.png";
[Ryy_est_true_mean, Ryy_est_samp_mean] = validate_periodogram_mean(hn, ryy_true, num_realizations, Ny, num_samp, true, savepath);

% estimate sample variance of the periodogram estimate
savepath = "../plots/prob1e/periodogram_var_comparison.png";
[Ryy_est_true_var, Ryy_est_samp_var] = validate_periodogram_var(hn, Ryy_est_true_mean, Ryy_true, num_realizations, Ny, num_samp, true, savepath);


% validate the statistical properties of the periodogram estimate
% generate several realization of the periodogram estimate
num_realizations = 1000;
Ny = 4096*2+1;

num_samp_lst = transpose([10, 50:50:5000]);
len_lst = size(num_samp_lst, 1);

periodogram_mean_arr = zeros(len_lst, num_realizations);
periodogram_var_arr = zeros(len_lst, num_realizations);

disp("Computing mean and variance across 1000 realizations...");
for i = 1:num_realizations
    [xn, yn] = generate_random_process(hn, Ny, false);
    for j = 1:len_lst
        num_samp = num_samp_lst(j,1);
        y_samp = yn((Ny+1)/2 - num_samp/2 : (Ny+1)/2 + num_samp/2 - 1, 1);
        assert(size(y_samp,1) == num_samp);

        Ryy_est = estimate_periodogram(y_samp, false, "");
        periodogram_mean_arr(j, i) = mean(Ryy_est);
        periodogram_var_arr(j, i) = var(Ryy_est);
    end
end

fig = figure;
plot(num_samp_lst, mean(periodogram_mean_arr, 2), LineWidth=2); hold on;
plot(num_samp_lst, mean(real(Ryy_true)) + num_samp_lst * 0, LineWidth=2); hold off;
xlabel("number of samples N", FontSize=16);
ylabel("Periodogram mean", FontSize=16);
title("Periodogram Mean over 1000 realizations", FontSize=16);
grid on;
legend("sample mean", "true mean");
saveas(fig, "../plots/prob1e/periodogram_sample_mean.png");
close;

fig = figure;
plot(num_samp_lst, mean(periodogram_var_arr, 2), LineWidth=2); hold on;
plot(num_samp_lst, var(real(Ryy_true)) + num_samp_lst * 0, LineWidth=2); hold off;
xlabel("number of samples N", FontSize=16);
ylabel("Periodogram variance", FontSize=16);
title("Periodogram Variance over 1000 realizations", FontSize=16);
grid on;
legend("sample variance", "true variance");
saveas(fig, "../plots/prob1e/periodogram_sample_var.png");
close;