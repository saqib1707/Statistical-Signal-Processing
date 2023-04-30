clc; clear;

[xn, Fs] = audioread("../data/x.wav");
[v2n, ~] = audioread("../data/v2.wav");  % highly correlated with v1[n]; no speech
Nx = size(xn,1);
Nv2 = size(v2n,1);

% sound(xn,Fs);
% sound(v2n,Fs);

fig = figure;
plot(0:Nx-1, xn);
xlabel("sample index [n]", FontSize=14);
ylabel("amplitude", FontSize=14);
title("Noisy speech signal x[n]", FontSize=16);
grid on;
saveas(fig, "../plots/noisyAudio.png");
close;

fig = figure;
plot(0:Nv2-1, v2n);
xlabel("sample index [n]", FontSize=14);
ylabel("amplitude", FontSize=14);
title("Noise signal v_2[n]", FontSize=16);
grid on;
saveas(fig, "../plots/v2Noise.png");
close;

% rxnxn = estimate_autocorr(xn,true,"../plots/noisySpeechAutocorr.png");
rv2v2 = estimate_autocorr(v2n,true,"../plots/v2NoiseAutocorr.png");
rxnv2 = estimate_crosscorr(xn,v2n,true,"../plots/noisySpeechv2NoiseCrossCorr.png");
% rv2xn = estimate_crosscorr(v2n,xn,true,"../plots/noisySpeechv2NoiseCrossCorr.png");

% rv2v2 = xcorr(v2n,v2n);
% rxnv2 = xcorr(xn,v2n);

% rxnv2_cpy = zeros(2*Nx-1,1);
% rxnv2_cpy(1:Nx-1) = flip(rv2xn(Nx+1:2*Nx-1));
% rxnv2_cpy(Nx,1) = rv2xn(Nx,1);
% rxnv2_cpy(Nx+1:2*Nx-1) = flip(rv2xn(1:Nx-1));
% assert( sum(abs(rxnv2 - conj(rxnv2_cpy)), 'all') == 0 );

filtOrd = [4,8,12,16,20];
for M = filtOrd

    RXX = compute_corrmat(rv2v2, M);   % (M,M)
    % RYX = compute_corrmat(rxnv2, M);
    RYX = transpose(rxnv2(Nx:Nx+M-1));  % (L,M)
    
    Co = (RYX * inv(RXX))';         % (M,1)
    disp("Filter coefficients M="+M);
    disp(Co');
    v1nhat1 = zeros(Nx,1);
    for i=1:Nx
        if i < M
            X = transpose([flip(v2n(1:i))', zeros(1,M-i)]);
        else
            X = flip(v2n(i-M+1:i));
        end
        v1nhat1(i) = Co' * X;
    end

    v1nhat2 = filter(Co,1,v2n);
    assert(sum(abs(v1nhat1 - v1nhat2), 'all') < 1e-10);
    %  v1nhat = rescale(v1nhat, -1, 1);
    
    fig = figure;
    plot(0:Nx-1, v1nhat1);
    xlabel("sample index [n]", FontSize=14);
    ylabel("amplitude", FontSize=14);
    title("Estimated noise v_1[n] (Filter order M="+M+")", FontSize=16);
    grid on;
    saveas(fig, "../plots/estv1n_M"+M+".png");
    close;
    
    snhat = (xn - v1nhat1);
%     snhat = rescale(snhat, -1, 1);
%     sound(snhat,Fs);
    audiowrite("../results/cleanAudio_M"+M+".wav", snhat, Fs);
    
    fig = figure;
    plot(0:Nx-1, snhat);
    xlabel("sample index [n]", FontSize=14);
    ylabel("amplitude", FontSize=14);
    title("Estimated clean audio s[n] (Filter order M="+M+")", FontSize=16);
    grid on;
    saveas(fig, "../plots/cleanAudio_M"+M+".png");
    close;
    
    % snhat2 = smoothdata(xn, "gaussian", 500);
    % sound(snhat2,Fs);
    % fig = figure;
    % plot(0:Nx-1, snhat2);
    % xlabel("samples [n]");
    % ylabel("amplitude");
    % grid on;
    % saveas(fig, "../plots/cleanSpeech2.png");
    % close;

    % compute SNR
    snrVal = snr(snhat, v1nhat1);
    disp("SNR based on estimated audio and estimated noise for M="+M+":"+snrVal);
    snrVal2 = 10 * log10(sum(snhat.^2) / sum(v1nhat1.^2));
    disp("SNR based on estimated audio and estimated noise for M="+M+":"+snrVal2);

    % compute variance
    varVal = var(snhat);
    disp("Variance of the estimated audio for M="+M+":"+varVal);
end