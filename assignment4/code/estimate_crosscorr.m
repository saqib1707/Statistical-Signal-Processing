function [rxy] = estimate_crosscorr(xn, yn, tosave, savepath)
    % estimate autocorrelation sequence using N samples
    Nx = size(xn,1);
    Ny = size(yn,1);
    assert(Nx == Ny);
    N = Nx;

    rxy = zeros(2*N-1,1);    % 2(N-1)+1
    for m = 0:N-1
        rxy(N+m,1) = sum(xn(1+m:N,1) .* conj(yn(1:N-m,1)));
    end
    for m = -(N-1):-1
        rxy(N+m,1) = sum(xn(1:N-abs(m), 1) .* conj(yn(abs(m)+1:N, 1)));
    end
    rxy = rxy / N;

    if tosave == true
        fig = figure;
        plot(-(N-1):(N-1), rxy, LineWidth=1);
    %     plot(-(N-1):(N-1), rxx2, LineWidth=2); hold on;
        xlabel("lag [m]", FontSize=14);
        ylabel("Cross-correlation r_{xy}[m]", FontSize=14);
        title("Estimated Cross-correlation ("+N+" samples)", FontSize=14);
        grid on;
        saveas(fig, savepath);
        close;
    end

end