clc;
clear;
close all

%% Plot the total intensity over time of phage-treated immunocompetent mice

% load mice image data
dir_path = '../data/';
files = dir([dir_path '*.mat']);

for f = 1:length(files)
    load([dir_path files(f).name])
    disp(files(f).name)
end

% only work with I+/P+ data
IposPpos_data = [data_IposPpos_2015_02_24_female data_IposPpos_2015_02_24_male data_RagPpos_2015_02_24 data_RagPpos_2011_09_20_mice123 data_RagPpos_2011_09_20_mice456];
time_labels = {'H2'; 'H4'; 'H6'; 'H8'; 'H24'; 'H48'; 'H72'};

% Calculate total intensity signal for top and bottom compartments
[intensity_top, intensity_bottom] = calculate_total_intensity(IposPpos_data, time_labels);

IposPpos_pixint_top = intensity_top;
IposPpos_pixint_bottom = intensity_bottom;


indx_alive = find(~isnan(IposPpos_pixint_top(:,end)));
live_top = IposPpos_pixint_top(indx_alive, :);
live_bottom = IposPpos_pixint_bottom(indx_alive, :);
times = [2 4 6 8 24 48 72];
[clearance_top, clearance_bottom, indx_significant] = find_cleartimes_test(live_top, live_bottom, 3, times);

%Plot total intensity over time for WT/P+ mice
figure(1);
xticks_lab  = string([ 2 4 6 8 24 48 72]) + 'H';
yl = [0 30];
title_lab = {'Top'};
t = tiledlayout(2, 1, 'TileSpacing','Compact');
nexttile
plot_fun_signal_intensity(live_top(indx_significant,:), xticks_lab, yl, title_lab, 'Total intensity', 0, 60, 50)
hold on
plot([1 7], [3 3], '--k', 'LineWidth', 1.5, 'DisplayName', 'Intensity threshold')
hold off
nexttile
title_lab = {'Bottom'};
plot_fun_signal_intensity(live_bottom(indx_significant,:), xticks_lab, yl, title_lab, 'Total intensity', 0, 60, 50)
hold on
plot([1 7], [3 3], '--k', 'LineWidth', 1.5, 'DisplayName', 'Intensity threshold')
hold off
set(gcf, 'position', [1001         746         842         591])
% file_name = '../figures/Fig8a_total_intensity_allmice.eps';
% exportgraphics(gcf, file_name)

%% save total intensity signal of the two compartments
immuno_data_top = IposPpos_pixint_top;
immuno_data_bottom = IposPpos_pixint_bottom;
% save immuno_data_top.mat immuno_data_top
% save immuno_data_bottom.mat immuno_data_bottom