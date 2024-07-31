function M = gps_position(data)
    %% interpolation of GPS location
    dataGPS = data.NMEA.string;
    
    M = [];
    for i=1:size(dataGPS,2)
        data_i = dataGPS{i};
        parsed_data_i = parseNMEA(data_i);
        if parsed_data_i.type == "GNGLL"
            n = size(M,1);
            M(n+1, 1) = parsed_data_i.lat;
            M(n+1, 2) = parsed_data_i.lat_hem;
            M(n+1, 3) = parsed_data_i.lon;
            M(n+1, 4) = parsed_data_i.lon_hem;
        end
    end
    % warning, cols 2 and 4 are actually chars!
    % char(87)=='W'
    % char(83)=='S'

end

% fprintf("Creating GPS file\n");
% gps_locations = interp1(1:size(M,1), M, linspace(1,size(M,1),size(H,2)));
% fileID = fopen([nb{1}(1:end-4) '.gpx'], 'W');
% fprintf(fileID, '<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\" ?>\n');
% fprintf(fileID, '<gpx version=\"1.1\" creator=\"CustomGPX\" xmlns=\"http://www.topografix.com/GPX/1/1\">\n');
% for i=1:size(H,2)
%     fprintf(fileID, '<wpt lat=\"-%.4f\" lon=\"-%.4f\"><name>x</name></wpt>\n', gps_locations(i,1), gps_locations(i,3));
% end
% fprintf(fileID, '</gpx>\n');
% fclose(fileID);
% fprintf("Done\n");
% save(sprintf("%s_f_%d.mat", nb{1}(1:end-4), fband), "H", "gps_locations", "depths", "PulseDuration", "PulseForm", "TransmitPower");
