%% Calculate infection clearance time from data of phage-treated immunocompetent mice

clc
clear
close all

%% find clearance time for immunocompetent mice that survived

% load data for I+/P+ mice
load('../data/immuno_data_top.mat')
load('../data/immuno_data_bottom.mat')

indx_alive = find(~isnan(immuno_data_bottom(:,end)));

live_top = immuno_data_top(indx_alive, :);
live_bottom = immuno_data_bottom(indx_alive, :);

times = [2 4 6 8 24 48 72];
% function to calculate infection clearance time using total intensity
% signal
[clearance_top, clearance_bottom, indx_significant] = find_cleartimes_test(live_top, live_bottom, 3, times);

% compared infection clerance time of bottom vs top compartment
[p,h,stats] = signrank(clearance_bottom, clearance_top, 'tail', 'left')


figure(1);
plot_columns_clearance_IposCase(clearance_top, clearance_bottom, [],'Clearance time (h)', p)
% file_name = '../figures/Fig8b_time_clearance_allmice.eps';
% exportgraphics(gcf, file_name)
