clc
clear
close all

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');  

list1=ls ('../BroadBandProcessing/BURRIED_OBJECTS_2/D20240420-T15*.raw', '-1');
list2=ls ('../BroadBandProcessing/BURRIED_OBJECTS_2/D20240420-T16*.raw', '-1');
list = [list1 list2];
list = splitlines(list);
[longList,~]=size(list);

% RAW DATA
figure('Position', [675 573 570 2*413]);
tiledlayout(3,1);
for archiv= [2 5 6]
    
    filename=strtrim(list(archiv,:));

    [header,data]=readEK80(filename);
    cell = list(archiv);
    path = cell{1};
    name = path(end-20:end-4);
    load(['MAT/' name '_f_1.mat']);


    sig_lf = data.pings(1).comp_sig_1;
    [ntimelf,npingslf] = size(sig_lf);
    
    nexttile;
    imagesc(1:npingslf, (1:ntimelf)*params.SampleInterval * 1500/2, 20*log10(abs(sig_lf)));
    xlabel("Pings");
    ylabel("Depth (in m)");
    colorbar;
    title(sprintf("Low frequencies, Slope=%.2f, PulseDuration=%dμs", params.Slope, 1e6*params.PulseDuration));
    if archiv == 2
        xlim([100 500]);
    end
    if archiv == 5
        xlim([125 525]);
    end
    if archiv == 6
        xlim([200 600]);
    end
end
saveas(gcf, 'FIG/lf_raw.eps', 'epsc');
drawnow;

% LOW FREQUENCIES
figure('Position', [675 573 570 2*413]);
tiledlayout(3,1);
    
for archiv= [2 5 6]%1:longList-1

    filename=strtrim(list(archiv,:));

    [header,data]=readEK80(filename);
    cell = list(archiv);
    path = cell{1};
    name = path(end-20:end-4);
    load(['MAT/' name '_f_1.mat']);

    sig_lf = data.pings(1).comp_sig_1;
    [ntimelf,npingslf] = size(sig_lf);
    
    
    nexttile;
    imagesc(1:size(H,2), f/1e3, abs(H));
    xlabel("Pings");
    ylabel("Frequency (in kHz)");
    colorbar;
    if archiv == 2
        xlim([100 500]);
    end
    if archiv == 5
        xlim([125 525]);
    end
    if archiv == 6
        xlim([200 600]);
    end
    title(sprintf("Low frequencies, Slope=%.2f, PulseDuration=%dμs", params.Slope, 1e6*params.PulseDuration));
end
saveas(gcf, 'FIG/lf_256.eps', 'epsc');



% HIGH FREQUENCIES
figure('Position', [675 573 570 2*413]);
tiledlayout(3,1);
    
for archiv= [2 5 6]%1:longList-1

    filename=strtrim(list(archiv,:));

    [header,data]=readEK80(filename);
    cell = list(archiv);
    path = cell{1};
    name = path(end-20:end-4);
    load(['MAT/' name '_f_2.mat']);
    sig_hf = data.pings(2).comp_sig_1;
    [ntimehf,npingshf] = size(sig_hf);

    nexttile;
    imagesc(1:size(H,2), f/1e3, abs(H));
    xlabel("Pings");
    ylabel("Frequency (in kHz)");
    colorbar;
    if archiv == 2
        xlim([100 500]);
    end
    if archiv == 5
        xlim([125 525]);
    end
    if archiv == 6
        xlim([200 600]);
    end
    title(sprintf("High frequencies, Slope=%.2f, PulseDuration=%dμs", data.params(2).Slope, 1e6*data.params(2).PulseDuration));
end

saveas(gcf, 'FIG/hf_256.eps', 'epsc');
