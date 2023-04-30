function [Ryy_est_true_mean, Ryy_est_samp_mean] = validate_periodogram_mean(hn, ryy_true, num_realizations, Ny, num_samp, tosave, savepath)
    % estimate sample periodogram mean over many realizations
    Ryy_est_samp_mean = zeros(2 * num_samp - 1, 1);
    
    for i = 1:num_realizations
        [xn, yn] = generate_random_process(hn, Ny, false);
    
        % y_samp = yn(1:num_samp, 1);
        y_samp = yn((Ny+1)/2 - num_samp/2 : (Ny+1)/2 + num_samp/2 - 1, 1);
        assert(size(y_samp,1) == num_samp);
    
        Ryy_est = estimate_periodogram(y_samp, false, "");
        Ryy_est_samp_mean = Ryy_est_samp_mean + Ryy_est;
    end
    
    Ryy_est_samp_mean = Ryy_est_samp_mean / num_realizations;
    
    % estimate actual periodogram mean using the formula in slides
    N_ryy_true = size(ryy_true, 1);
    m = -(N_ryy_true - 1)/2 : (N_ryy_true - 1)/2;
    
    N = 2 * num_samp - 1;
    Ryy_est_true_mean = zeros(N, 1);
    for n = 0 : N-1
        window = transpose(1 - abs(m) / ((N_ryy_true + 1)/2));
        exp_term = transpose(exp(-1j * 2 * pi * n * m / N));
        Ryy_est_true_mean(n+1, 1) = sum(window .* ryy_true .* exp_term);
    end

    if tosave == true
        fig = figure;
        plot((2*pi/N)*(0:N-1), Ryy_est_samp_mean, LineWidth=1); hold on;
        plot((2*pi/N)*(0:N-1), real(Ryy_est_true_mean), LineWidth=1); hold off;
        xlabel("\omega", FontSize=16);
        ylabel("Mean E(R_{yy}(e^{jw}))", FontSize=16);
        title("Comparison of sample and true mean of periodogram estimate", FontSize=16);
        grid on;
        legend("sample mean", "true mean");
        saveas(fig, savepath);
        close;
    end
end