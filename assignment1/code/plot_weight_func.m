function [] = plot_weight_func(weight_func, win_len, win_shift, tosave, prob, save_file_name)
    fig1 = figure;
    plot(0:size(weight_func, 1)-1, weight_func, 'LineWidth', 2.0);
    grid on;
    title("Weight Function for Hamming Win (M = "+win_len+", R = "+win_shift+")", 'FontSize', 18);
    xlabel("n", 'FontSize', 18);
    ylabel('$\tilde{w}[n]$', 'interpreter', 'latex', 'FontSize', 20, 'fontweight', 'bold');
    
    if tosave == true
        saveas(fig1, "../plots/prob"+prob+"/"+save_file_name);
        close;
    end

    fig2 = figure;
    plot(0:size(weight_func, 1)-1, weight_func, 'LineWidth', 2.0); hold on;
    plot(0:size(weight_func, 1)-1, hamming(win_len), 'LineWidth', 2.0); hold off;
    grid on;
    title("Weight Function and Hamming Win (M = "+win_len+", R = "+win_shift+")", 'FontSize', 18);
    xlabel("n", 'FontSize', 18);
    ylabel('signal intensity', 'interpreter', 'latex', 'FontSize', 20, 'fontweight', 'bold');
    legend("weight func", "hamming win")
    
    if tosave == true
        saveas(fig2, "../plots/prob"+prob+"/"+"weight_func_vs_hamming_M256_R128.png");
        close;
    end
end