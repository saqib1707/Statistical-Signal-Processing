clc; clear;

audio_filepath = "../data/exhaust.wav";
[audio_data, samp_freq] = audioread(audio_filepath);
num_samp = size(audio_data,1);
% sound(audio_data, samp_freq);

fig = figure;
plot((0:num_samp-1)/samp_freq, audio_data, LineWidth=1);
xlabel("time [s]", FontSize=16);
ylabel("amplitude", FontSize=16);
title("Exhaust Audio Signal", FontSize=16);
grid on;
saveas(fig, "../plots/prob2a/exhaust_audio_signal.png");
close;


% estimate autocorrelation of the audio signal
disp("Computing autocorrelation of audio signal...");
savepath = "../plots/prob2a/exhaust_autocorr.png";
ryy_audio = estimate_autocorr(audio_data, true, savepath);
disp("Estimated autocorrelation of audio signal");
 
% estimate periodogram of the audio signal
disp("Computing periodogram of audio signal...");
savepath = "../plots/prob2a/exhaust_periodogram.png";
Ryy_audio = estimate_periodogram(audio_data, true, savepath);
disp("Estimated periodogram of audio signal");
