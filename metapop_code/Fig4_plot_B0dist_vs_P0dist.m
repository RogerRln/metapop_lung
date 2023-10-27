%% Code to simulate the effects of different forms of allocating phage dose and bacterial inoculum in the network
% Inoculum: Phage-sensitive bacteria (BP)
% Phage added two hours after infection

% To simulate from scratch the different distributions of phage dose and
% bacterial inoculum uncomment lines 35 to 248

clc
clear;
close all; 

lung_mass = 0.135; % lung mass in grams

% Immunocompetence and phage therapy parameters
B = 1e6; % CFU, total bacteria
P = 1e7; % PFU, total phage
I = (2.7e6*lung_mass); % cells, initial number of neutrophils (from Roach et al., 2017)

% maximum number of neutrophils, we assume a large pool of neutrophils is available (10X Nlung)
nlung = 3.24e+06;
max_neutrophil_num = nlung*10; 

%% Uniform distribution of bacterial inoculum and phage dose

% bacterial inoculum and phage dose distribution
b_dist = 1:15; % Uniform distribution of bacterial inoculum
p_dist = 1:15; % Uniform distribution of phage dose

simu_time = 100; % simulation time in hours

% Simulate the metapopulation model
[time, res, p] = simu_metapop(B, P, I, max_neutrophil_num, b_dist, p_dist, simu_time);
data_uniBac_uniPhage.res = res;
data_uniBac_uniPhage.time = time;
%save('../data/data_uniBac_uniPhage.mat', 'data_uniBac_uniPhage')

% %% Uniform distribution of bacterial inoculum and top distribution of phage dose
% 
% % bacterial inoculum and phage dose distribution
% b_dist = 1:15; % Uniform distribution of bacterial inoculum
% p_dist = 1:3; % adding phage in the first 3 nodes
% 
% simu_time = 100; % simulation time in hours
% 
% % Simulate the metapopulation model
% [time, res, p] = simu_metapop(B, P, I, max_neutrophil_num, b_dist, p_dist, simu_time);
% data_uniBac_topPhage.res = res;
% data_uniBac_topPhage.time = time;
% save('../data/data_uniBac_topPhage.mat', 'data_uniBac_topPhage')
% 
% %% Uniform distribution of bacterial inoculum and bottom distribution of phage dose
% 
% % bacterial inoculum and phage dose distribution
% b_dist = 1:15; % Uniform distribution of bacterial inoculum
% p_dist = 4:15; % adding phage in the last 12 nodes
% 
% simu_time = 100; % simulation time in hours
% 
% % Simulate the metapopulation model
% [time, res, p] = simu_metapop(B, P, I, max_neutrophil_num, b_dist, p_dist, simu_time);
% data_uniBac_bottomPhage.res = res;
% data_uniBac_bottomPhage.time = time;
% save('../data/data_uniBac_bottomPhage.mat', 'data_uniBac_bottomPhage')
% 
% %% Uniform distribution of bacterial inoculum and adding phage to node 1
% 
% % bacterial inoculum and phage dose distribution
% b_dist = 1:15; % Uniform distribution of bacterial inoculum
% p_dist = 1; % adding phage in node 1
% 
% simu_time = 100; % simulation time in hours
% 
% % Simulate the metapopulation model
% [time, res, p] = simu_metapop(B, P, I, max_neutrophil_num, b_dist, p_dist, simu_time);
% data_uniBac_node1Phage.res = res;
% data_uniBac_node1Phage.time = time;
% save('../data/data_uniBac_node1Phage.mat', 'data_uniBac_node1Phage')
% 
% 
% %% Top distribution of bacterial inoculum and uniform distribution of phage dose
% 
% % bacterial inoculum and phage dose distribution
% b_dist = 1:3; % adding bacteria in the first 3 nodes
% p_dist = 1:15; % uniform distribution of phage
% 
% simu_time = 100; % simulation time in hours
% 
% % Simulate the metapopulation model
% [time, res, p] = simu_metapop(B, P, I, max_neutrophil_num, b_dist, p_dist, simu_time);
% data_topBac_uniPhage.res = res;
% data_topBac_uniPhage.time = time;
% save('../data/data_topBac_uniPhage.mat', 'data_topBac_uniPhage')
% 
% %% Top distribution of bacterial inoculum and phage dose
% 
% % bacterial inoculum and phage dose distribution
% b_dist = 1:3; %  adding bacteria in the first 3 nodes
% p_dist = 1:3; %  adding phage in the first 3 nodes
% 
% simu_time = 100; % simulation time in hours
% 
% % Simulate the metapopulation model
% [time, res, p] = simu_metapop(B, P, I, max_neutrophil_num, b_dist, p_dist, simu_time);
% data_topBac_topPhage.res = res;
% data_topBac_topPhage.time = time;
% save('../data/data_topBac_topPhage.mat', 'data_topBac_topPhage')
% 
% %% Top distribution of bacterial inoculum and bottom distribution of phage dose
% 
% % bacterial inoculum and phage dose distribution
% b_dist = 1:3;  %  adding bacteria in the first 3 nodes
% p_dist = 4:15; %  adding bacteria in the last 12 nodes
% 
% simu_time = 100; % simulation time in hours
% 
% % Simulate the metapopulation model
% [time, res, p] = simu_metapop(B, P, I, max_neutrophil_num, b_dist, p_dist, simu_time);
% data_topBac_bottomPhage.res = res;
% data_topBac_bottomPhage.time = time;
% save('../data/data_topBac_bottomPhage.mat', 'data_topBac_bottomPhage')
% 
% 
% %% Top distribution of bacterial inoculum and adding phage in node 1
% 
% % bacterial inoculum and phage dose distribution
% b_dist = 1:3;  %  adding bacteria in the first 3 nodes
% p_dist = 1; %  adding bacteria in node 1
% 
% simu_time = 100; % simulation time in hours
% 
% % Simulate the metapopulation model
% [time, res, p] = simu_metapop(B, P, I, max_neutrophil_num, b_dist, p_dist, simu_time);
% data_topBac_node1Phage.res = res;
% data_topBac_node1Phage.time = time;
% save('../data/data_topBac_node1Phage.mat', 'data_topBac_node1Phage')
% 
% %% Bottom distribution of bacterial inoculum and uniform distribution of phage dose
% 
% % bacterial inoculum and phage dose distribution
% b_dist = 4:15;  %  adding bacteria in the last 12 nodes
% p_dist = 1:15; %  uniform distribution of phage dose
% 
% simu_time = 100; % simulation time in hours
% 
% % Simulate the metapopulation model
% [time, res, p] = simu_metapop(B, P, I, max_neutrophil_num, b_dist, p_dist, simu_time);
% data_bottomBac_uniPhage.res = res;
% data_bottomBac_uniPhage.time = time;
% save('../data/data_bottomBac_uniPhage.mat', 'data_bottomBac_uniPhage')
% 
% %% Bottom distribution of bacterial inoculum and top distribution of phage dose
% 
% % bacterial inoculum and phage dose distribution
% b_dist = 4:15;  % adding bacteria in the last 12 nodes
% p_dist = 1:3; %   adding phage in the first 3 nodes
% 
% simu_time = 100; % simulation time in hours
% 
% % Simulate the metapopulation model
% [time, res, p] = simu_metapop(B, P, I, max_neutrophil_num, b_dist, p_dist, simu_time);
% data_bottomBac_topPhage.res = res;
% data_bottomBac_topPhage.time = time;
% save('../data/data_bottomBac_topPhage.mat', 'data_bottomBac_topPhage')
% 
% %% Bottom distribution of bacterial inoculum and phage dose
% 
% % bacterial inoculum and phage dose distribution
% b_dist = 4:15;  %  adding bacteria in the last 12 nodes
% p_dist = 4:15; %  adding phage in the last 12 nodes
% 
% simu_time = 100; % simulation time in hours
% 
% % Simulate the metapopulation model
% [time, res, p] = simu_metapop(B, P, I, max_neutrophil_num, b_dist, p_dist, simu_time);
% data_bottomBac_bottomPhage.res = res;
% data_bottomBac_bottomPhage.time = time;
% save('../data/data_bottomBac_bottomPhage.mat', 'data_bottomBac_bottomPhage')
% 
% %% Bottom distribution of bacterial inoculum and adding phage in node 1
% 
% % bacterial inoculum and phage dose distribution
% b_dist = 4:15;  % adding bacteria in the last 12 nodes
% p_dist = 1; %  adding phage in node 1
% 
% simu_time = 100; % simulation time in hours
% 
% % Simulate the metapopulation model
% [time, res, p] = simu_metapop(B, P, I, max_neutrophil_num, b_dist, p_dist, simu_time);
% data_bottomBac_node1Phage.res = res;
% data_bottomBac_node1Phage.time = time;
% save('../data/data_bottomBac_node1Phage.mat', 'data_bottomBac_node1Phage')
% 
% 
% %% Adding bacterial inoculum in node 1 and uniform distribution of phage dose
% 
% % bacterial inoculum and phage dose distribution
% b_dist = 1;  % adding bacteria in node 1
% p_dist = 1:15; %  uniform distribution of phage dose
% 
% simu_time = 100; % simulation time in hours
% 
% % Simulate the metapopulation model
% [time, res, p] = simu_metapop(B, P, I, max_neutrophil_num, b_dist, p_dist, simu_time);
% data_node1Bac_uniPhage.res = res;
% data_node1Bac_uniPhage.time = time;
% save('../data/data_node1Bac_uniPhage.mat', 'data_node1Bac_uniPhage')
% 
% %% Adding bacterial inoculum in node 1 and top distribution of phage dose
% 
% % bacterial inoculum and phage dose distribution
% b_dist = 1;  % adding bacteria in node 1
% p_dist = 1:3; %  adding phage in the first 3 nodes
% 
% simu_time = 100; % simulation time in hours
% 
% % Simulate the metapopulation model
% [time, res, p] = simu_metapop(B, P, I, max_neutrophil_num, b_dist, p_dist, simu_time);
% data_node1Bac_topPhage.res = res;
% data_node1Bac_topPhage.time = time;
% save('../data/data_node1Bac_topPhage.mat', 'data_node1Bac_topPhage')
% 
% %% Adding bacterial inoculum in node 1 and bottom distribution of phage dose
% 
% % bacterial inoculum and phage dose distribution
% b_dist = 1;  % adding bacteria in node 1
% p_dist = 4:15; %  adding phage in the last 12 nodes
% 
% simu_time = 100; % simulation time in hours
% 
% % Simulate the metapopulation model
% [time, res, p] = simu_metapop(B, P, I, max_neutrophil_num, b_dist, p_dist, simu_time);
% data_node1Bac_bottomPhage.res = res;
% data_node1Bac_bottomPhage.time = time;
% save('../data/data_node1Bac_bottomPhage.mat', 'data_node1Bac_bottomPhage')
% 
% %% Adding bacterial inoculum and phage dose in node 1
% 
% % bacterial inoculum and phage dose distribution
% b_dist = 1;  % adding bacteria in node 1
% p_dist = 1; %  adding phage in node 1
% 
% simu_time = 100; % simulation time in hours
% 
% % Simulate the metapopulation model
% [time, res, p] = simu_metapop(B, P, I, max_neutrophil_num, b_dist, p_dist, simu_time);
% data_node1Bac_node1Phage.res = res;
% data_node1Bac_node1Phage.time = time;
% save('../data/data_node1Bac_node1Phage.mat', 'data_node1Bac_node1Phage')


%% Load data of the different phage dose and bacterial inoculum distributions

load '../data/data_uniBac_uniPhage.mat'
load '../data/data_uniBac_topPhage.mat'
load '../data/data_uniBac_bottomPhage.mat'
load '../data/data_uniBac_node1Phage.mat'

load '../data/data_topBac_uniPhage.mat'
load '../data/data_topBac_topPhage.mat'
load '../data/data_topBac_bottomPhage.mat'
load '../data/data_topBac_node1Phage.mat'

load '../data/data_bottomBac_uniPhage.mat'
load '../data/data_bottomBac_topPhage.mat'
load '../data/data_bottomBac_bottomPhage.mat'
load '../data/data_bottomBac_node1Phage.mat'

load '../data/data_node1Bac_uniPhage.mat'
load '../data/data_node1Bac_topPhage.mat'
load '../data/data_node1Bac_bottomPhage.mat'
load '../data/data_node1Bac_node1Phage.mat'


% load data_uniBac_uniPhage.mat
% load data_uniBac_topPhage.mat
% load data_uniBac_bottomPhage.mat
% load data_uniBac_node1Phage.mat
% 
% load data_topBac_uniPhage.mat
% load data_topBac_topPhage.mat
% load data_topBac_bottomPhage.mat
% load data_topBac_node1Phage.mat
% 
% load data_bottomBac_uniPhage.mat
% load data_bottomBac_topPhage.mat
% load data_bottomBac_bottomPhage.mat
% load data_bottomBac_node1Phage.mat
% 
% load data_node1Bac_uniPhage.mat
% data_node1Bac_uniPhage = data_node1Bac_uniformPhage;
% load data_node1Bac_topPhage.mat
% load data_node1Bac_bottomPhage.mat
% load data_node1Bac_node1Phage.mat

all_data = [data_uniBac_uniPhage data_topBac_uniPhage data_bottomBac_uniPhage data_node1Bac_uniPhage...
    data_uniBac_topPhage data_topBac_topPhage data_bottomBac_topPhage data_node1Bac_topPhage...
    data_uniBac_bottomPhage data_topBac_bottomPhage data_bottomBac_bottomPhage data_node1Bac_bottomPhage...
    data_uniBac_node1Phage data_topBac_node1Phage data_bottomBac_node1Phage data_node1Bac_node1Phage];

%% Plot the bacterial dynamics across the metapopulation network

figure(1)
t = tiledlayout(4,4, 'TileSpacing', 'Compact');
when_cbar = 4:4:16;
when_xlabel = 13:16;
when_ylabel = [1 5 9 13];
when_title = 1:4;
Bacinocu_titles = {'Uniform distribution of $B_0$', 'Top distribution of $B_0$', 'Bottom distribution of $B_0$', 'Inoculating $B_0$ at Node 1'};
Phageinocu_titles = {'Uniform dist. of $P_0$', 'Top dist. of $P_0$', 'Bottom dist. of $P_0$', 'Inoculating $P_0$ at Node 1'};
count_v = 1;
for i = 1:size(all_data,2)
    nexttile(i)
    
    if ismember(i, when_title)
        plot_meatpopDyn_InocuDist(all_data(i).res, all_data(i).time, p, Bacinocu_titles{i})
    else
        plot_meatpopDyn_InocuDist(all_data(i).res, all_data(i).time, p, '')
    end
    if ismember(i, when_cbar)
        cbar = colorbar;
        caxis([0 10]);
        ylabel(cbar, 'log$_{10}\,$ Total bacteria (CFU/ml)', 'interpreter', 'latex', 'fontsize', 15)
    end
    if ismember(i, when_xlabel)
        xlabel('Time (h)', 'interpreter', 'latex', 'fontsize', 17)
    end
    if ismember(i, when_ylabel)
        ylabel('Node index, g', 'interpreter', 'latex', 'fontsize', 17)
        txt = text(-14, 7.5, Phageinocu_titles{count_v},  'HorizontalAlignment', 'center',...
            'rotation', 90,  'Fontweight', 'bold', 'FontSize', 20 , 'Interpreter', 'latex');
        count_v = count_v + 1;
    end
    
end
set(gcf, 'position', [665         276        1251        1045])
% filename = '/Users/rrodriguez77/Dropbox (GaTech)/Phage-Immune_host-pathogen project/draft_metapopulation/version2_Sep132022/figures/B0Dist_vs_P0Dist/bac_density_time.png';
% exportgraphics(gcf, filename, 'Resolution', 300)

%% Plot the phage population dynamics across the metapopulation network

figure(2)
t = tiledlayout(4,4, 'TileSpacing', 'Compact');
when_cbar = 4:4:16;
when_xlabel = 13:16;
when_ylabel = [1 5 9 13];
when_title = 1:4;
Bacinocu_titles = {'Uniform distribution of $B_0$', 'Top distribution of $B_0$', 'Bottom distribution of $B_0$', 'Inoculating $B_0$ at Node 1'};
Phageinocu_titles = {'Uniform dist. of $P_0$', 'Top dist. of $P_0$', 'Bottom dist. of $P_0$', 'Inoculating $P_0$ at Node 1'};
count_v = 1;
for i = 1:size(all_data,2)
    nexttile(i)
    
    if ismember(i, when_title)
        plot_phageDyn_InocuDist(all_data(i).res, all_data(i).time, p, Bacinocu_titles{i})
    else
        plot_phageDyn_InocuDist(all_data(i).res, all_data(i).time, p, '')
    end
    if ismember(i, when_cbar)
        cbar = colorbar;
        caxis([0 12]);
        ylabel(cbar, 'log$_{10}\,\,$ Phage (PFU/ml)', 'interpreter', 'latex', 'fontsize', 17)
    end
    if ismember(i, when_xlabel)
        xlabel('Time (h)', 'interpreter', 'latex', 'fontsize', 17)
    end
    if ismember(i, when_ylabel)
        ylabel('Node index, g', 'interpreter', 'latex', 'fontsize', 17)
        txt = text(-10, 7.5, Phageinocu_titles{count_v},  'HorizontalAlignment', 'center',...
            'rotation', 90,  'Fontweight', 'bold', 'FontSize', 20 , 'Interpreter', 'latex');
        count_v = count_v + 1;
    end
    
end
set(gcf, 'position', [665         276        1251        1045])
% filename = '/Users/rrodriguez77/Dropbox (GaTech)/Phage-Immune_host-pathogen project/draft_metapopulation/version2_Sep132022/figures/B0Dist_vs_P0Dist/phage_density_time.eps';
% exportgraphics(gcf, filename)


%% Plot the time to infection resolution given variations in the phage and bacteria inoculation site


data_uniBac = [data_uniBac_uniPhage data_uniBac_topPhage data_uniBac_bottomPhage data_uniBac_node1Phage];
data_topBac =  [data_topBac_uniPhage data_topBac_topPhage data_topBac_bottomPhage data_topBac_node1Phage];
data_bottomBac = [data_bottomBac_uniPhage data_bottomBac_topPhage data_bottomBac_bottomPhage data_bottomBac_node1Phage];
data_node1Bac = [data_node1Bac_uniPhage data_node1Bac_topPhage data_node1Bac_bottomPhage data_node1Bac_node1Phage];

node_indx = 1:15;
cols_lines = brewermap(size(data_uniBac,2), 'Dark2');
markers_lines = {'o', 'square', 'diamon', 'x'};
figure(3)
t = tiledlayout(2,2);
for i = 1:size(data_uniBac,2)
    
    t_clear = time_clearance_thresh(data_uniBac(i).res, data_uniBac(i).time, p, 1./p.branch_volume);
    nexttile(1)
    %plot(node_indx, t_clear, 'color', cols_lines(i,:), 'LineWidth', 2.5) % plot color lines
    plot(node_indx, t_clear,  markers_lines{i}, 'color', 'k', 'LineWidth', 1.5, 'MarkerSize', 10) %plot with markers instead on lines
    hold on
    
end
xlabel('Node index, g', 'Interpreter', 'latex')
ylabel('Clearance time (h)', 'Interpreter', 'latex')
xticks(1:2:15)
xlim([1 15])
ylim([28 38])
title('Uniform distribution of $B_0$', 'Interpreter', 'latex', 'fontsize', 20)
legend('Uniform distribution of $P_0$', 'Top distribution of $P_0$', 'Bottom distribution of $P_0$', 'Inoculating $P_0$ at Node 1', 'Interpreter','latex', 'fontsize', 17)
legend boxoff
set(gca, 'fontsize', 17, 'LineWidth', 1.5)

for i = 1:size(data_topBac, 2)
    
    t_clear = time_clearance_thresh(data_topBac(i).res, data_topBac(i).time, p, 1./p.branch_volume);
    nexttile(2)
    %plot(node_indx, t_clear, 'color', cols_lines(i,:), 'LineWidth', 2.5)
    plot(node_indx, t_clear,  markers_lines{i}, 'color', 'k', 'LineWidth', 1.5, 'MarkerSize', 10) %plot with markers instead on lines
    hold on
    
end
xlabel('Node index, g', 'Interpreter', 'latex')
ylabel('Clearance time (h)', 'Interpreter', 'latex')
xticks(1:2:15)
xlim([1 15])
ylim([28 38])
title('Top distribution of $B_0$', 'Interpreter', 'latex', 'fontsize', 20)
legend('Uniform distribution of $P_0$', 'Top distribution of $P_0$', 'Bottom distribution of $P_0$', 'Inoculating $P_0$ at Node 1', 'Interpreter','latex', 'fontsize', 17)
legend boxoff
set(gca, 'fontsize', 18, 'LineWidth', 1.5)

for i = 1:size(data_bottomBac, 2)
    
    t_clear = time_clearance_thresh(data_bottomBac(i).res, data_bottomBac(i).time, p, 1./p.branch_volume);
    nexttile(3)
    %plot(node_indx, t_clear, 'color', cols_lines(i,:), 'LineWidth', 2.5)
    plot(node_indx, t_clear,  markers_lines{i}, 'color', 'k', 'LineWidth', 1.5, 'MarkerSize', 10) %plot with markers instead on lines
    hold on
    
end
xlabel('Node index, g', 'Interpreter', 'latex')
ylabel('Clearance time (h)', 'Interpreter', 'latex')
xticks(1:2:15)
xlim([1 15])
ylim([28 38])
title('Bottom distribution of $B_0$', 'Interpreter', 'latex', 'fontsize', 20)
legend('Uniform distribution of $P_0$', 'Top distribution of $P_0$', 'Bottom distribution of $P_0$', 'Inoculating $P_0$ at Node 1', 'Interpreter','latex', 'fontsize', 17)
legend boxoff
set(gca, 'fontsize', 18, 'LineWidth', 1.5)

for i = 1:size(data_node1Bac, 2)
    
    t_clear = time_clearance_thresh(data_node1Bac(i).res, data_node1Bac(i).time, p, 1./p.branch_volume);
    nexttile(4)
    %plot(node_indx, t_clear, 'color', cols_lines(i,:), 'LineWidth', 2.5)
    plot(node_indx, t_clear,  markers_lines{i}, 'color', 'k', 'LineWidth', 1.5, 'MarkerSize', 10) %plot with markers instead on lines
    hold on
    
end
xlabel('Node index, g', 'Interpreter', 'latex')
ylabel('Clearance time (h)', 'Interpreter', 'latex')
xticks(1:2:15)
xlim([1 15])
ylim([28 38])
title('Inoculating $B_0$ at Node 1', 'Interpreter', 'latex', 'fontsize', 20)
legend('Uniform distribution of $P_0$', 'Top distribution of $P_0$', 'Bottom distribution of $P_0$', 'Inoculating $P_0$ at Node 1', 'Interpreter','latex', 'fontsize', 17, 'location', 'southeast')
legend boxoff
set(gca, 'fontsize', 18, 'LineWidth', 1.5)
set(gcf, 'position', [1001         876         939         713], 'renderer', 'painters')

% filename = '/Users/rrodriguez77/Dropbox (GaTech)/Phage-Immune_host-pathogen project/draft_metapopulation/version2_Sep132022/figures/B0Dist_vs_P0Dist/clearance_time_InocuSites.png';
% exportgraphics(gcf, filename, 'Resolution', 300)
