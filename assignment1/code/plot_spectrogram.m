function [] = plot_spectrogram(spec_single_sided, freq_axis, time_axis, win_len, win_shift, tosave, prob, save_file_name)
    
    fig = figure;
    imagesc(time_axis, freq_axis, spec_single_sided);
    set(gca, 'ydir', 'normal');   % flip the Y Axis so lower frequencies are at the bottom
    
    title("Spectrogram (M = "+win_len+", R = "+win_shift+")", 'FontSize', 18);
    xlabel("Time [s]", 'FontSize', 16);
    ylabel("Frequency [KHz]", 'FontSize', 16);

    colormap(jet(256));
    c = colorbar;
    c.Label.String = 'Power/Freq [dB/Hz]';
    c.Label.FontSize = 16;

%     spec_max = max(spec_single_sided(:));
%     spec_min = min(spec_single_sided(:));
%     disp(spec_min);
%     disp(spec_max);
%     set(c, 'ylim', [spec_min spec_max]);

    if tosave == true
        saveas(fig, "../plots/prob"+prob+"/"+save_file_name);
        close;
    end

end