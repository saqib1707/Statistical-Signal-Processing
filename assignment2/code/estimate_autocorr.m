function [rxx1] = estimate_autocorr(xn, tosave, savepath)
    % estimate autocorrelation sequence using N samples
    N = size(xn, 1);
    rxx1 = zeros(2 * N - 1, 1);    % 2(N-1)+1
    
    for m = 0:N-1
        rxx1(N+m,1) = sum(xn(1+m:N,1) .* conj(xn(1:N-m,1)));
    end
    rxx1(1:N-1,1) = flip(rxx1(N+1:2*N-1, 1));
    rxx1 = rxx1 / N;

%     rxx2 = zeros(2 * N - 1, 1);    % 2(N-1)+1
%     for m = 0:N-1
%         rxx2(N+m,1) = sum(xn(1+m:N,1) .* conj(xn(1:N-m,1)));
%     end
%     for m = -(N-1):-1
%         rxx2(N+m,1) = sum(xn(1:N-abs(m), 1) .* conj(xn(abs(m)+1:N, 1)));
%     end
%     rxx2 = rxx2 / N;
%     disp(sum(abs(rxx1 - rxx2)));

    if tosave == true
        fig = figure;
        plot(-(N-1):(N-1), rxx1, LineWidth=1);
    %     plot(-(N-1):(N-1), rxx2, LineWidth=2); hold on;
        xlabel("index [n]", FontSize=16);
        ylabel("Autocorrelation r_{yy}[m]", FontSize=16);
        title("Estimated Autocorrelation ("+N+" samples)", FontSize=16);
        grid on;
        saveas(fig, savepath);
        close;
    end
end