function [t,u] = precomputation_pulse(data, freq_range, expected_duration)
    params =data.params(freq_range);
    PulseDuration = params.PulseDuration;
    SampleInterval = params.SampleInterval;
      
    pings = data.pings(freq_range);
    signals = mean(cat(3, pings.comp_sig_1, pings.comp_sig_2, pings.comp_sig_3), 3, 'omitnan');
    signal = signals(:,1); % we use the first signal

    signal0 = signal(1:floor(3*PulseDuration/SampleInterval));
    energy0 = signal0.*conj(signal0);
    centroid = floor(((1:length(energy0))*energy0) / sum(energy0));
    index_chirp = floor(centroid - expected_duration/SampleInterval/2) : floor(centroid + expected_duration/SampleInterval/2);
    positive_index_chirp = index_chirp(index_chirp>0);
    u = signal(positive_index_chirp);
    t = linspace(0,length(u) * SampleInterval, length(u));
end

% function [t, u] = precomputation_pulse(data, fband)
% 
%     params =data.params(fband);
%     PulseDuration = params.PulseDuration;
%     Slope = params.Slope;
% 
%     config = data.config(fband);
%     fmin = config.FrequencyMinimum;
%     fmax = config.FrequencyMaximum;
% 
%     stages = data.filter_coeff(fband).stages;
% 
%     fs = 1.5e6;
%     t = 0:1/fs:PulseDuration;
%     s = analytic_signal(fs,PulseDuration,fmin,fmax,Slope);
%     % % TODO simon check if this changes something?
%     % s = [zeros(1, 100) s zeros(1,100)];
%     % t0 = linspace(0,PulseDuration + 2000/fs, length(s));
%     [u,t] = filtering(s,t0,stages);
%     % normalization
%     u = u/max(abs(u));
%     u = u.';
% end
% 
% function s = analytic_signal(fs, PulseDuration,fmin,fmax,Slope)
%     time_chirp = 0:1/fs:PulseDuration;
%     s = chirp(time_chirp,fmin,PulseDuration, fmax);
%     weight = hann(2*round(Slope * fs * PulseDuration));
%     for i=1:round(Slope * fs * PulseDuration)
%         s(i) = weight(i) * s(i);
%         s(end+1-i) = weight(end+1-i) * s(end+1-i);
%     end
% end
% 
% function [ss, tt] = filtering(s,t,stages)
%     ss = s;
%     tt = t;
%     f_s_dec = zeros(length(stages),1);
%     for stage = 1:length(stages)
%         F = stages(stage);
%         filt = (F.Coefficients(1:2:end)+1i*F.Coefficients(2:2:end));
%         D = F.DecimationFactor;
%         if stage == 1
%             f_s_dec(stage) = 1/(t(2)-t(1))/D;
%         else
%             f_s_dec(stage) = f_s_dec(stage-1)/D;
%         end
%         % normalization of the signal
%         ss = ss / max(abs(ss), [],  'omitnan');
%         % filter
%         ss = conv(ss, filt,'same');
%         % down sampling
%         ss = downsample(ss, D);
%         tt = downsample(tt, D);
%     end
% end
% 
