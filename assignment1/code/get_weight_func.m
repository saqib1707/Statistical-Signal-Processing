function [weight_func] = get_weight_func(win_len, win_shift, win_type)
    if win_type == "hamming"
        win_samp = hamming(win_len);
    end

    win_shift_right_R = [zeros(win_shift, 1); win_samp(1:win_len-win_shift, 1)];
    win_shift_left_R = [win_samp(win_shift+1:win_len, 1); zeros(win_shift, 1)];
    
    if win_shift == win_len / 2
        weight_func = win_samp + win_shift_right_R + win_shift_left_R;
    end
end