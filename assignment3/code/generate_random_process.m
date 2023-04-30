function [xn, yn] = generate_random_process(hn, Ny, tosave)
    Nh = size(hn,1);
    Nx = Ny - Nh + 1;
    xn = randn(Nx, 1);

    if tosave == true
        fig = figure;
        plot(0:Nx-1, xn, LineWidth=1);
        xlabel("time index [n]", FontSize=14);
        ylabel("amplitude x[n]", FontSize=14);
        title("Gaussian White Noise Realization ~ N(0,1)", FontSize=14);
        grid on;
        saveas(fig, "../plots/gaussian_white_input_"+Ny+".png");
        close;
    end
    
    % estimate the random process y[n] by convolving x[n] and h[n]
    yn = conv(xn, hn);
    assert(size(yn,1) == Ny);

    if tosave == true
        fig = figure;
        plot(0:Ny-1, yn, LineWidth=1);
        xlabel("time index [n]", FontSize=14);
        ylabel("amplitude y[n]", FontSize=14);
        title("Output random process realization y[n]", FontSize=14);
        grid on;
        saveas(fig, "../plots/out_random_process_"+Ny+".png");
        close;
    end
end