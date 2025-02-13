function [nmea,nmea_type]=parseNMEA(nmea_string)
nmea_string(isspace(nmea_string))=' ';
idx = strfind(nmea_string, ',');

if isempty(idx)
    nmea=struct();
    nmea_type='invalid';
    return;
end

type = upper(nmea_string(2:idx(1) - 1));
nmeadata = nmea_string(idx(1) + 1:end);
%  remove checksum - add trailing comma for output lacking last field
idx = strfind(nmeadata, '*');

if (isempty(idx))
    nmeadata = [nmeadata ','];
else
    nmeadata = [nmeadata(1:idx - 1) ','];
end

nmeadata = strrep(nmeadata, ',,', ',0');

%     '$PRDID,0.14,0.97,83.17  '
%     '$INGST,065803.372,,0.8,0.5,5.6,0.6,0.7,0.8  '


switch type(3:end)
    case 'SHR'

        if contains(nmeadata,'ATT')
            %'$PASHR,ATT,348466.00,147.26310,-0.66647,0.28198,0.0011,0.0037,0 ,';
            format = '%s %2.0f %2d %f %f %f %f %f %f %f %f %d %d';
            out = textscan(nmeadata, format, 1, 'delimiter', ',');
            nmea = struct('type', type, ...
                'time', [out{2} out{3} out{4}], ...
                'heading', double(out{5}), ...
                'heading_type', 'T', ...
                'roll', out{6}, ...
                'pitch', out{7}, ...
                'heave', out{8});
            nmea_type='attitude';
        else
            %'$PASHR,065803.372,83.17,T,0.97,0.14,-0.36,0.021,0.021,0.015,2,1 ';
            format = '%2.0f %2d %f %f %c %f %f %f %f %f %f %d %d';
            out = textscan(nmeadata, format, 1, 'delimiter', ',');
            nmea = struct('type', type, ...
                'time', [out{1} out{2} out{3}], ...
                'heading', double(out{4}), ...
                'heading_type', out{5}, ...
                'roll', out{6}, ...
                'pitch', out{7}, ...
                'heave', out{8});
            nmea_type='attitude';
        end

    case 'HDT'
        %HEHDT,048.15,T*17
        format = '%f';
        out = textscan(nmeadata, format, 1, 'delimiter', ',');
        nmea = struct('type', type, ...
            'heading', out{1});
        nmea_type='heading';
    case 'GGA'
        %$GPGGA,120926,6913.4620,S,17811.5690,W,1,09,00.9,26.2,M,-59.6,M,,*6B
        %$INGGA,120926.645,6913.46117,S,17811.55897,W,2,09,1.2,-6.59,M,,,28,0000*0D
        %$GPGGA,010901.00,3850.99301424,S,17828.18369031,E,5,07,1.8,19.801,M,21.853,M,9.0,
        %$GPGGA,         ,5420.0279,N,01154.7914,E,,,,,,,,,

        switch nmeadata(1)
            case ',' %when gps does not have time in GPGGA
                format = '%2d %f %c %3d %f %c %d %d %f %f %c %f %c %f %d';
                out = textscan(nmeadata, format, 1, 'delimiter', ',');
                nmea = struct('type', type, ...
                    'time', [nan nan nan], ...
                    'lat', double(out{1}) + out{2} / 60, ...
                    'lat_hem', out{3}, ...
                    'lon', double(out{4}) + out{5} / 60, ...
                    'lon_hem', out{6});
            otherwise
                format = '%2d %2d %f %2d %f %c %3d %f %c %d %d %f %f %c %f %c %f %d';
                out = textscan(nmeadata, format, 1, 'delimiter', ',');

                %         nmea = struct('type', type, ...
                %             'time', [out{1} out{2} out{3}], ...
                %             'lat', double(out{4}) + out{5} / 60, ...
                %             'lat_hem', out{6}, ...
                %             'lon', double(out{7}) + out{8} / 60, ...
                %             'lon_hem', out{9}, ...
                %             'fix', out{10}, ...
                %             'nsat', out{11}, ...
                %             'precision', out{12}, ...
                %             'msl_alt', out{13}, ...
                %             'msl_unit', out{14}, ...
                %             'geoidal_alt', out{15}, ...
                %             'geoidal_unit', out{16}, ...
                %             'dif_age', out{17}, ...
                %             'dif_sta', out{18} ...
                %             );
                nmea = struct('type', type, ...
                    'time', [out{1} out{2} out{3}], ...
                    'lat', double(out{4}) + out{5} / 60, ...
                    'lat_hem', out{6}, ...
                    'lon', double(out{7}) + out{8} / 60, ...
                    'lon_hem', out{9});
        end
        nmea_type='gps';

    case 'RMC'

        %$GPRMC,022345.976,A,4118.1066,S,17448.3002,E,0.33,24.92,260116,,,A*44
        format = '%2.0f %2d %f %s %2d %f %c %3d %f %c %f %d %2d %2d %2d ';
        out = textscan(nmeadata, format, 1, 'delimiter', ',');
        nmea = struct('type', type, ...
            'time', [out{1} out{2} out{3}], ...
            'lat', double(out{5}) + out{6} / 60, ...
            'lat_hem', out{7}, ...
            'lon', double(out{8}) + out{9} / 60, ...
            'lon_hem', out{10},...
            'sog',out{11});
        nmea_type='gps';

    case 'GLL'
        %  parse the rest of the nmea text
        format = '%2d %f %c %3d %f %c %2d %2d %f %c';
        out = textscan(nmeadata, format, 1, 'delimiter', ',');

        %  convert status to GGA fix
        fix = strcmpi('A', out{10});

        %  define GPS geographic position datagram
        nmea = struct('type', type, ...
            'lat', double(out{1}) + out{2} / 60, ...
            'lat_hem', out{3}, ...
            'lon', double(out{4}) + out{5} / 60, ...
            'lon_hem', out{6}, ...
            'time', [out{7} out{8} out{9}], ...
            'fix', fix ...
            );
        nmea_type='gps';
    case 'VTG'

        %  parse the rest of the nmea text
        format = '%f %c %f %c %f %c %f %c';
        out = textscan(nmeadata, format, 1, 'delimiter', ',');

        %  define Course Over Ground datagram
        nmea = struct('type', type, ...
            'true_cov', out{1}, ...
            'tcov_label', out{2}, ...
            'mag_cov', out{3}, ...
            'mcov_label', out{4}, ...
            'sog_knts', out{5}, ...
            'sogn_unit', out{6}, ...
            'sog_kph', out{7},...
            'sogk_unit', out{8} ...
            );
        nmea_type='speed';
    case 'VLW'
        %  parse the rest of the nmea text
        format = '%f %c %f %c';
        out = textscan(nmeadata, format, 1, 'delimiter', ',');


        nmea = struct('type', type, ...
            'total_cum_dist', out{1}, ...
            'tcd_unit', out{2}, ...
            'dist_since_reset', out{3}, ...
            'dsr_unit', out{4} ...
            );
        nmea_type='dist';
    case 'DFT'
        % $KMDFT,%.2f,M'
        format = '%f %c';
        out = textscan(nmeadata, format, 1, 'delimiter', ',');
        nmea = struct('type', type, ...
            'depth', out{1}, ...
            'unit', out{2}...
            );
        nmea_type='depth';
    otherwise
        %  unknown datagram type
        nmea = struct('type', type, 'string', nmeadata);
        nmea_type='unknown';

end

end