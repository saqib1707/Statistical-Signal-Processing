function [rxx_fft] = estimate_periodogram(xn, tosave, savepath)
    % estimate periodogram using N samples
    Nx = size(xn, 1);
    xn_fft = fft(xn, 2 * Nx - 1);
    rxx_fft = (1/Nx) * abs(xn_fft).^2;
    N_rxx = size(rxx_fft, 1);

    if tosave == true
        fig = figure;
        plot((2 * pi / N_rxx) * (0 : N_rxx - 1), rxx_fft, LineWidth=1);
        xlabel("\omega", FontSize=16);
        ylabel("Periodogram R(e^{jw})", FontSize=16);
        title("Estimated Periodogram ("+Nx+" samples)", FontSize=16);
        grid on;
        saveas(fig, savepath);
        close;
    end
end