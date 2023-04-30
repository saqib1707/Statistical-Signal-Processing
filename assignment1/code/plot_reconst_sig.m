function [] = plot_reconst_sig(x, x_reconst, Fs, win_len, win_shift, tosave, prob)
    % plot the original and reconstructed signal
    fig1 = figure;
    time_axis = (1:size(x,1)) / Fs;
    plot(time_axis, x); hold on;
    plot(time_axis, x_reconst); hold off;
    grid on;
    title("Original vs Reconstructed Signal", "FontSize", 18);
    xlabel("Time [s]", "FontSize", 16);
    ylabel("Signal Amplitude y[n]", "FontSize", 16);
    legend("x[n]", "x_r[n]", "FontSize", 14);
    
    if tosave == true
        saveas(fig1, "../plots/prob"+prob+"/origvsrecons_M"+win_len+"_R"+win_shift+".png");
    end

    fig2 = figure;
    time_axis = (1:size(x,1)) / Fs;
    
    subplot(1,2,1);
    plot(time_axis, x);
    ylim([-1.0 1.0]);
    grid on;
    title("Original Signal", "FontSize", 18);
    xlabel("Time [s]", "FontSize", 16);
    ylabel("Signal Amplitude y[n]", "FontSize", 16);

    subplot(1,2,2);
    plot(time_axis, x_reconst);
    ylim([-1.0 1.0]);
    grid on;
    title("Reconstructed Signal", "FontSize", 18);
    xlabel("Time [s]", "FontSize", 16);
    ylabel("Signal Amplitude y[n]", "FontSize", 16);

    pos = get(gcf, 'Position');
    set(gcf, 'Position', pos+[300 -50 300 -50]);

    if tosave == true
        saveas(fig2, "../plots/prob"+prob+"/orig_recons_M"+win_len+"_R"+win_shift+".png");
        close all;
    end
end