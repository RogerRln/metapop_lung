%% Code to simulate the effects of different mucin levels on the phage clearance time of a P.a. infection
% Inoculum: Phage-sensitive bacteria (BP)
% Phage added two hours after infection

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

% bacterial inoculum and phage dose distribution
b_dist = 1:15; % uniformly distribute bacteria among 15 nodes
p_dist = 1:15; % uniformly distribute phage among 15 nodes

simu_time = 100; % simulation time in hours
mucin_levels = [1 3 4]; % mucin concentrations in %

% save infection clearance time of network nodes for different mucin levels
num_nodes = 15;
t_clear_mat = zeros(num_nodes, numel(mucin_levels));

for i = 1:numel(mucin_levels)   
    mucin_conc = mucin_levels(i);
    % Simulate the metapopulation model
    [time, res, p] = simu_metapop_mucin(B, P, I, max_neutrophil_num, b_dist, p_dist, simu_time, mucin_conc);
    t_clear = time_clearance_thresh(res, time, p, 1./p.branch_volume);
    t_clear_mat(:, i) = t_clear;
end

%% Plot infection clearance times across the network

figure(1);
cols = brewermap(length(mucin_levels), 'Paired');
plot(t_clear_mat(:,1), 'o', 'color', 'k', 'LineWidth', 2, 'markersize', 10)
hold on
plot(t_clear_mat(:,2), 'square', 'color', 'k','LineWidth', 2, 'markersize', 10)
plot(t_clear_mat(:,3), 'diamond', 'color', 'k','LineWidth', 2, 'markersize', 10)
hold off
xlabel('Node index, g', 'Interpreter', 'latex')
xlim([1 15])
ylim([28 33])
xticks(1:15)
yticks(28:33)
ylabel('Clearance time (h)', 'Interpreter', 'latex')
legend(string(mucin_levels) + '% mucin', 'Interpreter', 'latex', 'Location', 'best', 'fontsize', 16)
legend box off
title('Infection clearance time', 'interpreter', 'latex')
set(gca, 'fontsize', 17, 'linewidth', 1.5);
set(gcf, 'renderer', 'painters')
% file_name = '/Users/rrodriguez77/Dropbox (GaTech)/Phage-Immune_host-pathogen project/draft_metapopulation/version2_Sep132022/figures/time_to_clearance/uniform_distribution/time_clearance_thresh_1CFU_1-4_mucin.png';
% exportgraphics(gcf,file_name, 'resolution', 300);