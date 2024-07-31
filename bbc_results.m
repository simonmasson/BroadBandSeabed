clc
clear
close all

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');  

list=ls ('../BroadBandProcessing/PAN_AZUCAR/*.raw', '-1');
list = splitlines(list);
[longList,~]=size(list);

% storing transfer function
H = [];
GPS = [];
lengths = [];
for archiv= 1:longList-1

    cell = list(archiv);
    path = cell{1};
    name = path(end-20:end-4);
    L = load(['MAT/' name '_f_1.mat']);
    H = [H L.H];
    GPS = [GPS L.GPS];
    lengths = [lengths size(H,2)];
end
f=L.f;

%% Write into a CSV file
ToCSV = [sum(abs(H).^2);GPS];
writematrix(ToCSV.', 'bbc.csv');
ToCSVFiltered = ToCSV(:, ToCSV(1,:)<1);
writematrix(ToCSVFiltered(:,1:50:end).', 'bbc_filt.csv');

%%

pos = [
    45.0455 65.7999;
    45.0469366666667 65.79679;
    45.046605 65.8126233333333;
    45.0483766666667 65.8123683333333;
    45.0553333333333 65.811415;
    45.058965 65.8116566666667;
    45.0663116666667 65.8128833333333;
    45.0770566666667 65.8138916666667;
    45.0815183333333 65.8097366666667;
    45.073981 65.809821;
    45.052145 65.798531;
    45.045094 65.798082;
    45.06177 65.801495;
    45.08125 65.798478;
    45.078352 65.804698;
    45.078234 65.81167;
    45.061054 65.807171;
    45.070411 65.801393;
    45.057015 65.806864;
    45.057652 65.797707
    ].';

Hlabel = zeros(128, size(pos,2));
index_label = [];
for i = 1:size(pos,2)
    dist = getDistanceFromLatLon(pos(:,i), GPS);
    [~,m] = min(abs(dist));
    Hlabel(:,i) = H(:,m);
    index_label = [index_label m];
end

figure;
plot(smooth(sum(abs(H).^2),100));
for l = lengths
    xline(l);
end


function d = getDistanceFromLatLon(p1, p2)
    dLat = p2(1,:) - p1(1);
    dLon = p2(2,:) - p2(2);
    a = sind(dLat/2) .* sind(dLat/2) + ...
        cosd(p1(1)) .* cosd(p2(1)) .* ...
        sind(dLon/2) .* sind(dLon/2);
    d = atan2(sqrt(a), sqrt(1 - a));
end