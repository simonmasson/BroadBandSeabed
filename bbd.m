clc
clear
close all

%% Loading the data
addpath('ESP3');

list1=ls ('./RAW/D20240420-T15*.raw', '-1');
list2=ls ('./RAW/D20240420-T16*.raw', '-1');
list = [list1 list2];
list = splitlines(list);
[longList,~]=size(list);

for archiv= 1:longList-1
    filename=strtrim(list(archiv,:));

    [header,data]=readEK80(filename);
    
    M = gps_position(data);
    
    %% Parsing the data
    for fband = 1%1:2

        Ndft = 128;
        ping_start = 1;
        ping_end = 0; % 0=end

        n_mean = 15;

        [h, f, space, time, crop_ruys, depths] = impulse_response(data, Ndft, fband, n_mean, ping_start, ping_end, false);
    end
    
    %%
    H = zeros(Ndft,length(h));
    for i = 1:length(h)
        H(:,i) = fft(h{i}, Ndft);
    end
    h2 = ifft(H, [], 1);
    %% crop_ruys in a matrix
    m = 0;
    for i = 1:length(crop_ruys)
        if size(crop_ruys{i},1) > m
            m = size(crop_ruys{i},1);
        end
    end
    MatRuys = zeros(length(crop_ruys), m);
    for i = 1:length(crop_ruys)
        for j = 1:size(crop_ruys{i},1)
            MatRuys(i,j) = crop_ruys{i}(j);
        end
    end

    GPS_lat = interp1(linspace(0,1,size(M,1)), M(:,1), linspace(0,1,size(H,2)));
    GPS_lon = interp1(linspace(0,1,size(M,1)), M(:,3), linspace(0,1,size(H,2)));
    GPS = [GPS_lat;GPS_lon];

    %% MAT EXPORT
    fprintf("Creating MAT file\n");
    params = data.params(fband);
    save("MAT/" + filename{1}(end-20:end-4) + "_f_" + fband + ".mat", "f", "H", "GPS", "params");
    fprintf("Done\n");

    % %% GPX EXPORT
    % fprintf("Creating GPS file\n");
    % fileID = fopen(['EXPORT/GPX/' filename{1}(end-20:end-4) '.gpx'], 'W');
    % fprintf(fileID, '<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\" ?>\n');
    % fprintf(fileID, '<gpx version=\"1.1\" creator=\"CustomGPX\" xmlns=\"http://www.topografix.com/GPX/1/1\">\n');
    % for i=1:size(GPS,2)
    %     fprintf(fileID, '<wpt lat=\"-%.8f\" lon=\"-%.8f\"><name>x</name></wpt>\n', GPS(1,i), GPS(2,i));
    % end
    % fprintf(fileID, '</gpx>\n');
    % fclose(fileID);
    % fprintf("Done\n");

    %% Plot
    figure;
    tiledlayout(2,1);
    
    % Cross correlation
    nexttile;
    imagesc(ping_start:ping_end, space, 20*log10(abs(MatRuys.')));
    colorbar;
    xlim([ping_start, size(H,2)]);
    xlabel("Ping number");
    ylabel("Depth (in m)");
    title("Matched filter");
    
    % Transfer function
    nexttile;
    imagesc(ping_start:ping_start + size(H,2), f/1e3, abs(H/Ndft));
    axis xy;
    ylabel("Frequency (in kHz)");
    xlabel("Ping number");
    colorbar;
    title("Impulse response in frequency domain");
    
    drawnow;
    saveas(gcf, "FIG/" + filename{1}(end-20:end-4) + "_f_" + fband + ".eps", 'epsc');

end
