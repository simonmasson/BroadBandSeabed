clc
clear
close all

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');  

%% Loading the data
addpath('ESP3');

filename = 'RAW/D20240420-T160412.raw';
[header,data]=readEK80(filename);

fband = 1;
fmin = data.params(fband).FrequencyStart;
fmax = data.params(fband).FrequencyEnd;

% params
params = data.params(fband);
SoundVelocity = params.SoundVelocity;
PulseDuration = params.PulseDuration;
SampleInterval = params.SampleInterval;
Slope = params.Slope;
SamplingRate = 1/SampleInterval;

%% Filters and decimation
stages = data.filter_coeff(fband).stages;
Fs = 1.5e6;
t0 = 0:1/Fs:PulseDuration;

% filter 1
F1 = stages(1);
filt1 = F1.Coefficients(1:2:end)+1i*F1.Coefficients(2:2:end);
D1 = F1.DecimationFactor;
% filter 2
F2 = stages(2);
filt2 = F2.Coefficients(1:2:end)+1i*F2.Coefficients(2:2:end);
D2 = F2.DecimationFactor;

% in frequency domain
Filt1 = fft(filt1);
freq_filt1 = linspace(-Fs/2,Fs/2, length(filt1));

Filt2 = fft(filt2);
freq_filt2 = linspace(-Fs/D1/2, Fs/D1/2, length(filt2));
while freq_filt2(end) < fmin
    freq_filt2 = freq_filt2 + Fs/D1;
end

figure;
tiledlayout(2,1);
nexttile;
hold on;
xline(fmin/1000, '--');
xline(fmax/1000, '--');
plot(freq_filt1/1000, abs(fftshift(Filt1)));
plot(freq_filt2/1000, abs(fftshift(Filt2)));
xlabel("Frequency (in kHz)");
ylabel("Modulus");
legend(["" "" "Filter 1" "Filter 2"]);
xlim([-20 100]);

nexttile;
hold on;
xline(fmin/1000, '--');
xline(fmax/1000, '--');
plot(freq_filt1/1000, angle(fftshift(Filt1)));
plot(freq_filt2/1000, angle(fftshift(Filt2)));
xlabel("Frequency (in kHz)");
ylabel("Phase");
legend(["" "" "Filter 1" "Filter 2"]);
xlim([-20 100]);
saveas(gcf, 'FIG/filters.eps', 'epsc');

%% Theoretic transmitted pulse
% Ideal chirp
s0 = chirp(t0,fmin,PulseDuration, fmax);
weight = hann(2*round(Slope * Fs * PulseDuration));
for i=1:round(Slope * Fs * PulseDuration)
    s0(i) = weight(i) * s0(i);
    s0(end+1-i) = weight(end+1-i) * s0(end+1-i);
end
% filter 1
s1 = conv(s0, filt1,'same');
% down sampling 1
s1 = downsample(s1, D1);
t1 = downsample(t0, D1);
% filter 2
s2 = conv(s1, filt2,'same');
% down sampling 2
s2 = downsample(s2, D2);
t2 = downsample(t1, D2);

figure;
tiledlayout(3,1);
nexttile;
hold on;
plot(1e3*t0, s0);
xlabel("Time (in ms)");
ylabel("Amplitude");
title("Analytic chirp");

nexttile;
hold on;
plot(1e3*t1, abs(s1));
plot(1e3*t1, real(s1));
plot(1e3*t1, imag(s1));
xlabel("Time (in ms)");
ylabel("Amplitude");
legend(["Modulus" "Real part" "Imaginary part"]);
title("Signal after one filter+downsampling step");

nexttile;
hold on;
plot(1e3*t2, abs(s2));
plot(1e3*t2, real(s2));
plot(1e3*t2, imag(s2));
xlabel("Time (in ms)");
ylabel("Amplitude");
legend(["Modulus" "Real part" "Imaginary part"]);
title("Output signal of the echosounder");
saveas(gcf, 'FIG/theo_u.eps', 'epsc');

%% Received signal
s = data.pings(fband).comp_sig_1(:,400);
u = s(11:49);
time = linspace(0,SampleInterval * length(s),length(s));

Ndft = 256;
f = linspace(0, SamplingRate, Ndft);
while f(end) < fmin 
    f = f + SamplingRate;
end
assert(f(1) < fmin);
assert(f(end) > fmax);

S = fft(s,Ndft);

figure;
tiledlayout(2,1);
nexttile;
plot(time*1e3, 20*log10(abs(s)));
xlabel("Time (in ms)");
ylabel("Power level (in dB)");
nexttile;
hold on;
xline(fmin/1e3,'--');
xline(fmax/1e3, '--');
plot(f/1e3, abs(S/Ndft));
xlabel("Frequency (in kHz)");
ylabel("Power (arbitrary units)");
saveas(gcf, 'FIG/rec_time_freq.eps', 'epsc');

%% Frequency components
window = 2 * SamplingRate * PulseDuration;
hann_window = hann(window);
overlap = window-1;

[ssf,ffp,ttp] = stft(s, SamplingRate, Window=hann_window, OverlapLength=overlap, FFTLength=Ndft, FrequencyRange='twosided');
ssf = ssf/Ndft; % we use [0 Fs] format here
figure;
tiledlayout(2,1);
nexttile;
imagesc(ttp*1e3, f/1e3,20*log10(abs(ssf)));
xlabel("Time (in ms)");
ylabel("Frequency (in kHz)");
title("Amplitude of the frequency components");
colorbar;
nexttile;
imagesc(ttp*1e3, f/1e3,unwrap(angle(ssf), [], 2));
xlabel("Time (in ms)");
ylabel("Frequency (in kHz)");
title("Phase of the frequency components");
colorbar;
saveas(gcf, 'FIG/rec_freq_compo.eps', 'epsc');
%% Phase correction
% U = fft(u,Ndft);
% green_functions = U .* exp(2*1i*pi*f.'*(ttp-ttp(1)).');
% phase_correction = conj(green_functions)./abs(green_functions);
% 
% figure;
% hold on;
% plot(ttp*1e3, unwrap(angle(ssf(100,:))));
% plot(ttp*1e3, unwrap(angle(green_functions(100,:)), [], 2));
% xlabel("Time (in ms)");
% ylabel("Unwrapped phase");
% legend(["Signal for f = 39kHz" "Theoretic phase of the green function"]);

%% Spreading correction
space = time*SoundVelocity/2;
time_stft = time(window/2:end-window/2);
space_stft = time_stft*SoundVelocity/2;
spreading_correction = space_stft.^1.5;

%% Absorption correction
gamma = precomputation_absorption(data, 20, f);
absorption_correction = exp(2 * gamma.' .* space_stft);

%% Corrected signal
corrected_ssf = ssf .* spreading_correction .* absorption_correction;

[corrected_s, ti] = istft(corrected_ssf,SamplingRate,'Window',hann_window,'OverlapLength',overlap,'FFTLength',Ndft, FrequencyRange='twosided');
corrected_s = corrected_s*Ndft;

figure('Position', [675 573 570 320]);
hold on;
plot(time*1e3, 20*log10(abs(s)));
plot(time*1e3, 20*log10(abs(corrected_s)));
xlabel("Time (in ms)");
ylabel("Power level (in dB)");
legend(["Raw signal" "Corrected signal"]);
saveas(gcf, 'FIG/ch3_4_comp_corr.eps', 'epsc');

%% Echo finding
crop_pulse = round(5 * 2 / SoundVelocity * SamplingRate);
y = s;
% correlations
ruy = xcorr(y,u);
ruy = ruy(length(y):end);
ruu = xcorr(u);
ruu = ruu(length(u):end);
% find peaks
ruy_except_pulse = ruy(crop_pulse:end);   
[max1,b] = max(abs(ruy_except_pulse));
peak1 = crop_pulse+b;
depth = space(peak1);
peak1_before = peak1;
while abs(ruy(peak1_before)) > max1/10e2 && peak1_before > 0
    % -60dB of energy, i.e. 1000 Pa
    peak1_before = peak1_before-1;
    disp(abs(ruy(peak1_before)));
end
peak1_after = peak1;
while abs(ruy(peak1_after)) > max1/10e2 && peak1_after < length(ruy)
    % -60dB of energy, i.e. 1000 Pa
    peak1_after = peak1_after+1;
end
crop_index = peak1_before:peak1_after;
peak1m60dB = peak1;
while abs(ruy(peak1m60dB)) > max1/10e3 && peak1m60dB < length(ruy)% -60dB of energy, i.e. 1000 Pa
    peak1m60dB = peak1m60dB+1;
end
if depth < 65
    [~,b] = max(abs(ruy(peak1+crop_pulse:end)));
    peak2 = peak1+crop_pulse+b;
end
% crop_index = peak1- PulseDuration*SamplingRate : peak1m60dB;
crop_ruy = ruy(crop_index);
crop_y = y(crop_index);
h = deconv(crop_ruy, ruu, Method='least-squares', RegularizationFactor=0);

figure;
hold on;
plot(time*1e3, 20*log10(abs(ruy)));
% E1
rectangle('Position', [time(peak1+PulseDuration*SamplingRate)*1e3 -50 1e3*(time(peak2-PulseDuration*SamplingRate)-time(peak1+PulseDuration*SamplingRate)) 110], 'LineStyle', '--');
text(30, 65, 'E1');
% For H
rectangle('Position', [1e3*time(peak1_before) -55 1e3*(time(peak1_after)-time(peak1_before)) 120], 'LineStyle', '--', 'LineWidth', 2);
text(20, 70, 'Here');
% E2
rectangle('Position', [1e3*time(peak2-PulseDuration*SamplingRate) -50 1e3*(time(peak2-PulseDuration*SamplingRate)-time(peak1)) 110], 'LineStyle', '--');
text(48, 65, 'E2');
xlabel("Time (in ms)");
ylabel("Cross-correlation energy (in dB)");
saveas(gcf,'FIG/E1E2.eps', 'epsc');

%% Transfer function
% Naive transfer function
H1 = fft(crop_y,Ndft)./fft(u,Ndft);
% Using regularization
h2 = deconv(crop_y, u, Method='least-squares');
h2_tik = deconv(crop_y, u, Method='least-squares', RegularizationFactor = 100);
% Using matched filter
RUU = fft(ruu, Ndft);
H3 = fft(crop_ruy,Ndft)./RUU;
% h4 = deconv(crop_ruy, ruu, Method='least-squares'); % same as H3

%% Plot
figure('Position', [675 573 570 320]);
hold on;
plot(f/1e3, abs(H1));
plot(f/1e3,abs(fft(h2, Ndft)));
plot(f/1e3,abs(fft(h2_tik, Ndft)));
plot(f/1e3, abs(H3), 'Linewidth', 2);
% plot(f/1e3, abs(fft(h4, Ndft)),'o');
xlabel("Frequency (in kHz)");
ylabel("Transfer function approximations");
legend(["Naive transfer function" "Least square regularization" "Tikhonov regularization" "Using correlations"]);
saveas(gcf,'FIG/h.eps', 'epsc');