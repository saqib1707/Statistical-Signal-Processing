clc; clear;

input_mu = 0;
input_var = 1;

% compute the impulse function h[n] from the transfer function H(z)
b = [1, -0.9, 0.81];     % numerator coeffs
a = [1, -2.76, 3.809, -2.654, 0.924];    % denominator coeffs

[hn, tn] = impz(b, a);
Nh = size(hn, 1);

% plot the impulse response
fig = figure;
plot(tn, hn, LineWidth=1);
xlabel("time [n]", FontSize=14);
ylabel("amplitude h[n]", FontSize=14);
title("Impulse Response h[n]", FontSize=14);
grid on;
saveas(fig, "../plots/impulse_response_impz.png");
close;

num_eval_pts = 2048;

% compute the true system PSD for the given H(z) = D(z) / A(z)
n = 0:(num_eval_pts - 1);
omega = 2 * pi * n / num_eval_pts;
Dz = 1 - 0.9 * exp(-1j * omega) + 0.81 * exp(-2j * omega);
Az = 1 - 2.76 * exp(-1j * omega) + 3.809 * exp(-2j * omega) - 2.654 * exp(-3j * omega) + 0.924 * exp(-4j * omega);

Hz = Dz ./ Az;
true_sys_psd = input_var * abs(Hz).^2;

num_realizations = 1000;
Ny = 8193;

% estimate autocorrelation sequence of y[n] for various sample lengths
num_samp_lst = [64, 256, 1024, 4096];
filter_order_lst = [2, 8, 14];

for num_samp = num_samp_lst
    for P = filter_order_lst
        disp("Case - N="+num_samp+", P="+P);
        true_psd_est_model = zeros(num_realizations, num_eval_pts);

        fig1 = figure;
        for itr = 1:num_realizations
            % generate samples of white gaussian noise random process x[n] ~ N(0,1)
            % Using generated x[n], generate a realization of the process y[n] = h[n] *
            % x[n] of length Ny 
            if itr == 1 && P == 2
                [xn, yn] = generate_random_process(hn, Ny, true);
            else
                [xn, yn] = generate_random_process(hn, Ny, false);
            end

            % sample y[n] from the center
            y_samp = yn((Ny+1)/2 - num_samp/2 : (Ny+1)/2 + num_samp/2 - 1, 1);
            assert(size(y_samp, 1) == num_samp);

            if itr == 1 && P == 2
                savepath = "../plots/est_autocorr_"+num_samp+".png";
                ryy_est = estimate_autocorr(y_samp, true, savepath);
            else
                ryy_est = estimate_autocorr(y_samp, false, "");
            end

            % create the autocorrelation toeplitz matrix
            autocorr_toeplitz_mat = compute_autocorr_toeplitz_mat(ryy_est, P);
    
            % 1st method
%             rhs_vec = zeros(P+1, 1);
%             rhs_vec(1, 1) = input_var;
%             ap_par_vec = autocorr_toeplitz_mat \ rhs_vec;
%             ap_par_vec = ap_par_vec / ap_par_vec(1,1);

            % 2nd method
            rhs_vec = -conj(ryy_est(num_samp+1:num_samp+P));
            ap_vec = autocorr_toeplitz_mat(2:P+1, 2:P+1) \ rhs_vec;
            sigma_sq = ryy_est(num_samp) + sum(ryy_est(num_samp+1:num_samp+P) .* ap_vec);
            ap_par_vec = zeros(P+1, 1);
            ap_par_vec(1,1) = 1;
            ap_par_vec(2:P+1, 1) = ap_vec;

            % compute true psd of this estimated model H(z) = 1 / A(z)
            Az_hat = compute_Az(ap_par_vec, P, num_eval_pts);
            Hz_hat = 1 ./ Az_hat;

            true_psd_est_model(itr, :) = sigma_sq * abs(Hz_hat).^2;

            plot((0 : num_eval_pts - 1) * 2 * pi / num_eval_pts, true_psd_est_model(itr, :), LineWidth=1); hold on;
        end

        hold off;
        xlabel("angular freq omega (rad)", FontSize=14);
        ylabel("R_{yy}(e^{jw})", FontSize=14);
        title("True PSD of Estimated Model (Ny="+num_samp+", P="+P+")", FontSize=14);
        grid on;
        saveas(fig1, "../plots/true_psd_est_model_Ny"+num_samp+"_P"+P+".png");
        close;

        % part (b) and (c)
        mean_est_psd = mean(true_psd_est_model, 1);
        median_est_psd = median(true_psd_est_model, 1);
        var_est_psd = var(true_psd_est_model, 1);
        
        % plot the true system PSD, mean and median of all realizations for
        % each case
        fig = figure;
        freq_axis = (0 : (num_eval_pts - 1)) * 2 * pi / num_eval_pts;
        plot(freq_axis, true_sys_psd, LineWidth=2, Color="r"); hold on;
        plot(freq_axis, mean_est_psd, LineWidth=2, Color="g"); hold on;
        plot(freq_axis, median_est_psd, LineWidth=2, Color="b"); hold off;
        xlabel("angular freq omega (rad)", FontSize=14);
        ylabel("R_{yy}(e^{jw})", FontSize=14);
        title("True system PSD, Mean, Median (Ny="+num_samp+", P="+P+")", FontSize=14);
        legend('True PSD', 'Mean PSD', 'Median PSD');
        grid on;
        saveas(fig, "../plots/true_sys_psd_mean_med_Ny"+num_samp+"_P"+P+".png");
        close;

        % plot the variance of all realizations for each case
        fig = figure;
        freq_axis = (0 : (num_eval_pts - 1)) * 2 * pi / num_eval_pts;
        plot(freq_axis, var_est_psd, LineWidth=2);
        xlabel("angular freq omega (rad)", FontSize=14);
        ylabel("Variance", FontSize=14);
        title("Variance of all realizations (Ny="+num_samp+", P="+P+")", FontSize=14);
        grid on;
        saveas(fig, "../plots/var_psd_Ny"+num_samp+"_P"+P+".png");
        close;
    end
end

close all;