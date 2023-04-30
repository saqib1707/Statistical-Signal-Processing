function [Ryy_true] = compute_true_periodogram(ryy, tosave, savepath)
    N_ryy = size(ryy, 1);
    Ryy_true = zeros(N_ryy, 1);

    m_arr = -(N_ryy - 1) / 2 : (N_ryy - 1) / 2;
    
    for n = 0 : N_ryy - 1
        exp_term = transpose(exp(-1j * 2 * pi * n * m_arr / N_ryy));
        Ryy_true(n+1, 1) = sum(ryy(m_arr + (N_ryy - 1)/2 + 1, 1) .* exp_term);
    end

%     Ryy_true_2 = zeros(N_ryy, 1);
%     for n = 0 : N_ryy - 1
%         for m = -(N_ryy - 1) / 2 : (N_ryy - 1) / 2
%             Ryy_true_2(n+1, 1) = Ryy_true_2(n+1, 1) + ryy(m + (N_ryy - 1)/2 + 1, 1) * exp(-1j * 2 * pi * n * m / N_ryy);
%         end
%     end
%     Ryy_diff = max(abs(real(Ryy_true) - real(Ryy_true_2)));
%     disp("Diff between custom FFT and system FFT:"+Ryy_diff);

%     Ryy_fft = fft(circshift(ryy, 550));
%     Ryy_diff = sum(abs(real(Ryy_true) - abs(Ryy_fft)));
%     disp("Diff between custom FFT and system FFT:"+Ryy_diff);

    if tosave == true
        fig = figure;
        plot((2 * pi / N_ryy) * (0 : N_ryy - 1), real(Ryy_true), LineWidth=1);
    %     plot((2 * pi / N_ryy) * (0 : N_ryy - 1), abs(Ryy_fft), LineWidth=1); hold off;
        xlabel("\omega", FontSize=16);
        ylabel("R_{yy}(e^{jw})", FontSize=16);
        title("True Power Spectral Density R_{yy}(e^{jw})", FontSize=16);
        grid on;
        saveas(fig, savepath);
        close;
    end
end