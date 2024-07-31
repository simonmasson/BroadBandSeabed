% function [h,f,space, time, ruys, depths] = impulse_response(data, Ndft, fband, n_mean, meters_before, meters_after, ping_start, ping_end)
function [h,f,space, time, crop_ruys, depths] = impulse_response(data, Ndft, fband, n_mean, ping_start, ping_end, FIXED_PEAK1_BEFORE)
    % params
    params = data.params(fband);
    SoundVelocity = params.SoundVelocity;
    PulseDuration = params.PulseDuration;
    SampleInterval = params.SampleInterval;
    SamplingRate = 1/SampleInterval;

    % Pulse computation
    [t_chirp, u] = precomputation_pulse(data, fband, 2*PulseDuration);

    r_uu = xcorr(u);
    r_uu = r_uu(length(u):end);

    fmin = data.config(fband).FrequencyMinimum;
    fmax = data.config(fband).FrequencyMaximum;

    % pings
    pings = data.pings(fband);
    signals = cat(3,pings.comp_sig_1,pings.comp_sig_2,pings.comp_sig_3,pings.comp_sig_4);
    if fband == 2
        n_sectors = 1;
    else
        n_sectors = 3;
    end

    f = linspace(0, SamplingRate, Ndft);
    while f(end) < fmin 
        f = f + SamplingRate;
    end
    assert(f(1) < fmin);
    assert(f(end) > fmax);

    window = 2 * SamplingRate * PulseDuration;
    hann_window = hann(window);
    overlap = window-1;

    R_UU = fft(r_uu, Ndft);

    %% time and space
    time = linspace(0,SampleInterval * size(signals,1), size(signals,1));
    space = time*SoundVelocity/2;
    time_stft = time(window/2:end-window/2);
    space_stft = time_stft*SoundVelocity/2;

    spherical_loss = space_stft.^1.5;
    depth_estimation = 20; % for absorption estimation
    gamma = precomputation_absorption(data, depth_estimation, f);
    absorption_loss = exp(2 * gamma.' .* space_stft);

    % iterating on all pings
    if ping_end == 0
        ping_end = size(signals, 2);
    end
    
    corrected_s = zeros(length(time), ping_end-ping_start, n_sectors);
    for ping = ping_start:ping_end
        fprintf("Correction ping %.0f%% \n", (ping-ping_start)/(ping_end-ping_start)*100);
        % corrections
        sigs = squeeze(signals(:,ping, :));
        sigs(isnan(sigs)) = 0;
        [ssf,ffp,ttp] = stft(sigs, SamplingRate, Window=hann_window, OverlapLength=overlap, FFTLength=Ndft, FrequencyRange='twosided');
        corrected_ssf = ssf .* spherical_loss .* absorption_loss;
        [corrected_ss, ti] = istft(corrected_ssf,SamplingRate,Window=hann_window,OverlapLength=overlap,FFTLength=Ndft,FrequencyRange='twosided');
        corrected_s(:, ping-ping_start+1,:) = corrected_ss(:,1:n_sectors);
    end

    ruys = zeros(length(time), ping_end-ping_start-n_mean, n_sectors);
    crop_ruys = {};
    h = {};
    depths = zeros(ping_end-ping_start-n_mean, 1);
    crop_pulse = round(5 * 2 / SoundVelocity * SamplingRate); % 2*length(u);
    for i = 1:ping_end-ping_start - n_mean
        fprintf("Transfer function ping %.0f%%\n", i/(ping_end-ping_start-n_mean)*100);
        % mean over sectors
        corrected_s_mean_sectors = mean(corrected_s, 3, 'omitnan');
        % mean over pings
        y = mean(corrected_s_mean_sectors(:,i:i+n_mean),2);
        % matched filter
        ruy = xcorr(y,u);
        ruy = ruy(length(y):end);

        ruys(:, i) = ruy;
        % find echo
        ruy_except_pulse = ruy(crop_pulse:end);   
        [max1,b] = max(abs(ruy_except_pulse));
        peak1 = crop_pulse+b;
        depths(i,1) = space(peak1);
        if FIXED_PEAK1_BEFORE
            peak1_before = peak1 - 2*PulseDuration*SamplingRate;
        else
            peak1_before = peak1;
            while abs(ruy(peak1_before)) > max1/10e2 && peak1_before > 1
                % -40dB of energy, i.e. 100 Pa
                peak1_before = peak1_before-1;
            end
        end
        peak1_after = peak1;
        while abs(ruy(peak1_after)) > max1/10e2 && peak1_after < length(ruy)
            % -40dB of energy, i.e. 100 Pa
            peak1_after = peak1_after+1;
        end
        crop_index = peak1_before:peak1_after;
        crop_ruys{i} = ruy(crop_index);
        h{i} = ifft(fft(crop_ruys{i},Ndft)./R_UU);
    end
end
