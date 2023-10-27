%% Code to simulate the metapopulation model under different immune and phage conditions
% Inoculum: Phage-sensitive bacteria (BP), uniformly distributed
% If phage therapy is used, we add phage 2 h after bacterial infection

clc
clear;
close all; 

lung_mass = 0.135; % lung mass in grams
simu_time = 100; % simulation time in hours
b_dist = 1:15; % uniformly distribute bacteria among 15 nodes
p_dist = 1:15; % uniformly distribute phage among 15 nodes

% maximum allowed number of lung neutrophils, we assume a large pool of neutrophils is available (10X Nlung)
nlung = 3.24e+06;
max_neutrophil_num = nlung*10; 

%% Immunocompetence and phage therapy parameters
B = 1e6; % CFU, total bacteria
P = 1e7; % PFU, total phage
I = (2.7e6*lung_mass); % cells, initial number of neutrophils (from Roach et al., 2017)
% Simulate the metapopulation model
[time, res, p] = simu_metapop(B, P, I, max_neutrophil_num, b_dist, p_dist, simu_time);
UniformInocu_Ipos_Ppos.res = res;
UniformInocu_Ipos_Ppos.time = time;

%% Immunocompetence without phage therapy parameters
B = 1e6; % CFU, total bacteria
P = 0; % PFU, total phage
I = (2.7e6*lung_mass); % cells, initial number of neutrophils (from Roach et al., 2017)

% Simulate the metapopulation model
[time, res, p] = simu_metapop(B, P, I, max_neutrophil_num, b_dist, p_dist, simu_time);
UniformInocu_Ipos_Pneg.res = res;
UniformInocu_Ipos_Pneg.time = time;

%% Neutropenia with phage therapy parameters
B = 1e6; % CFU, total bacteria
P = 1e7; % PFU, total phage
I = 0; % cells

% Simulate the metapopulation model
[time, res, p] = simu_metapop(B, P, I, max_neutrophil_num, b_dist, p_dist, simu_time);
UniformInocu_Ineg_Ppos.res = res;
UniformInocu_Ineg_Ppos.time = time;


%% Neutropenia without phage therapy parameters
B = 1e6; % CFU, total bacteria
P = 0; % PFU, total phage
I = 0; % cells

% Simulate the metapopulation model
[time, res, p] = simu_metapop(B, P, I, max_neutrophil_num, b_dist, p_dist, simu_time);
UniformInocu_Ineg_Pneg.res = res;
UniformInocu_Ineg_Pneg.time = time;


%% Plot the bacterial infection dynamics across the metapopulation network

% load '../data/UniformInocu_Ineg_Pneg.mat'
% load '../data/UniformInocu_Ineg_Ppos.mat'
% load '../data/UniformInocu_Ipos_Pneg.mat'
% load '../data/UniformInocu_Ipos_Ppos.mat'

color_bar_label = 'log$_{10}\,\,$ Total bacteria (CFU/ml)';
x_label = 'Time (h)';
y_label = 'Node index, g';
figure(1)
t = tiledlayout(2,2);
nexttile
plot_meatpopDyn(UniformInocu_Ineg_Pneg.res, UniformInocu_Ineg_Pneg.time, p, 'Immune-/Phage-', [], y_label, x_label)
nexttile
plot_meatpopDyn(UniformInocu_Ipos_Pneg.res, UniformInocu_Ipos_Pneg.time, p, 'Immune+/Phage-', color_bar_label, [], x_label)
nexttile
plot_meatpopDyn(UniformInocu_Ineg_Ppos.res, UniformInocu_Ineg_Ppos.time, p, 'Immune-/Phage+', [], y_label, x_label)
nexttile
plot_meatpopDyn(UniformInocu_Ipos_Ppos.res, UniformInocu_Ipos_Ppos.time, p, 'Immune+/Phage+', color_bar_label, [], x_label)
set(gcf, 'position', [1001         455        1008         882])
% filename = '/Users/rrodriguez77/Dropbox (GaTech)/Phage-Immune_host-pathogen project/draft_metapopulation/version2_Sep132022/figures/metapop_dynamics_ImmunePhage/uniform_distribution/metapop_dynamics';
% exportgraphics(gcf, [filename '.eps']);
% exportgraphics(gcf, [filename '.png']);

%% Plot the phage dynamics across the metapopulation network
color_bar_label = 'log$_{10}\,$ Phage (PFU/ml)';
x_label = 'Time (h)';
y_label = 'Node index';
figure(2)
t = tiledlayout(1,2, 'TileSpacing', 'Compact');
nexttile
plot_phageDyn(UniformInocu_Ineg_Ppos.res, UniformInocu_Ineg_Ppos.time, p, 'Immune-/Phage+',  [], y_label, x_label)
nexttile
plot_phageDyn(UniformInocu_Ipos_Ppos.res, UniformInocu_Ipos_Ppos.time, p, 'Immune+/Phage+',  color_bar_label, [], x_label)
set(gcf, 'position', [ 772         445        1067         433])
% filename = '/Users/rrodriguez77/Dropbox (GaTech)/Phage-Immune_host-pathogen project/draft_metapopulation/version2_Sep132022/figures/metapop_dynamics_ImmunePhage/uniform_distribution/metapop_Phage_dynamics';
% exportgraphics(gcf, [filename '.eps']);
% exportgraphics(gcf, [filename '.png']);


%% Plot the inset of the Immune+/Phage+ case
% we show times when infection clears from network nodes

figure(3)
color_bar_label = 'log$_{10}\,\,$ Total bacteria (CFU/ml)';
x_label = 'Time (h)';
y_label = 'Node index, g';
plot_meatpopDyn_zoomIposPpos(UniformInocu_Ipos_Ppos.res, UniformInocu_Ipos_Ppos.time, p, 'Immune+/Phage+', color_bar_label, y_label, x_label)
% file_name = '/Users/rrodriguez77/Dropbox (GaTech)/Phage-Immune_host-pathogen project/draft_metapopulation/version2_Sep132022/figures/metapop_dynamics_ImmunePhage/uniform_distribution/PNG_08092023/zoom_IposPpos.png';
% exportgraphics(gcf, file_name, 'Resolution', 300);