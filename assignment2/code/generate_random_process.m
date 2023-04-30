function [xn, yn] = generate_random_process(hn, Ny, tosave)
    Nh = size(hn,1);
    Nx = Ny - Nh + 1;
    xn = randn(Nx, 1);

    if tosave == true
        fig = figure;
        plot(0:Nx-1, xn, LineWidth=1);
        xlabel("index [n]", FontSize=16);
        ylabel("x[n]", FontSize=16);
        title("Gaussian White Noise Realization N(0,1)", FontSize=16);
        grid on;
        saveas(fig, "../plots/prob1b/gaussian_white_input.png");
        close;
    end
    
    % estimate the random process y[n] by convolving x[n] and h[n]
    yn = conv(xn, hn);
    assert(size(yn,1) == Ny);

    if tosave == true
        fig = figure;
        plot(0:Ny-1, yn, LineWidth=1);
        xlabel("index [n]", FontSize=16);
        ylabel("y[n]", FontSize=16);
        title("Output random process y[n]", FontSize=16);
        grid on;
        saveas(fig, "../plots/prob1b/out_random_process.png");
        close;
    end
end