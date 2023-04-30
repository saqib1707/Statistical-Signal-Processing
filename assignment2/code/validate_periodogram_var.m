function [Ryy_est_true_var, Ryy_est_samp_var] = validate_periodogram_var(hn, Ryy_est_true_mean, Ryy_true, num_realizations, Ny, num_samp, tosave, savepath)
    % estimate sample periodogram variance over many realizations
    N_Ryy_est = 2 * num_samp - 1;
    Ryy_est_samp_var = zeros(N_Ryy_est, 1);
    
    for i = 1:num_realizations
        [xn, yn] = generate_random_process(hn, Ny, false);
    
        % y_samp = yn(1:num_samp, 1);
        y_samp = yn((Ny+1)/2 - num_samp/2 : (Ny+1)/2 + num_samp/2 - 1, 1);
        assert(size(y_samp, 1) == num_samp);
    
        Ryy_est = estimate_periodogram(y_samp, false, "");
        Ryy_est_samp_var = Ryy_est_samp_var + (Ryy_est - Ryy_est_true_mean).^2;
    end
    
    Ryy_est_samp_var = Ryy_est_samp_var / (num_realizations-1);

    % estimate true variance of the periodogram estimate
    N_Ryy_true = size(Ryy_true, 1);
    omega_arr = transpose((2*pi / N_Ryy_true) * (0 : N_Ryy_true - 1));

    Ryy_est_true_var = (Ryy_true.^2) .* (1 + (sin(N_Ryy_true * omega_arr) ./ (N_Ryy_true * sin(omega_arr))).^2);
    disp(size(Ryy_est_true_var));
    disp(size(Ryy_est_samp_var));

    if tosave == true
        fig = figure;
        plot((2*pi / N_Ryy_est) * (0 : N_Ryy_est - 1), real(Ryy_est_samp_var), LineWidth=2); hold on;
        plot((2*pi / N_Ryy_true) * (0 : N_Ryy_true - 1), real(Ryy_est_true_var), LineWidth=1); hold off;
        xlabel("\omega", FontSize=16);
        ylabel("Variance Var(R_{yy}(e^{jw}))", FontSize=16);
        title("Comparison of sample and true variance of periodogram estimate", FontSize=16);
        grid on;
        legend("sample variance", "true variance");
        saveas(fig, savepath);
        close;
    end
end