function [rxx] = compute_true_autocorr(hn, tosave, savepath)
    Nh = size(hn, 1);

    rxx = zeros(2 * Nh - 1, 1);

    for m = -(Nh - 1) : (Nh - 1)
        for p = 0 : Nh - 1
            k = m + p;
            if k >= 0 && k <= Nh - 1
                rxx(Nh + m, 1) = rxx(Nh + m, 1) + hn(k + 1, 1) * conj(hn(p + 1, 1));
            end
        end
    end

    fig = figure;
    plot(-(Nh-1):(Nh-1), rxx, LineWidth=1);
    xlabel("index [n]", FontSize=16);
    ylabel("Autocorrelation r_{yy}[m]", FontSize=16);
    title("True Autocorrelation r_{yy}[m]", FontSize=16);
    grid on;

    if tosave == true
        saveas(fig, savepath);
        close;
    end
end