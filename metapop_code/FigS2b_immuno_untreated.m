%% Code to simulate the infection scenario where host is immunocompetent but not phage-treated
% Inoculum: Phage-sensitive bacteria (BP)


clc
clear;
close all; 

lung_mass = 0.135; % lung mass in grams

% Immunocompetence with out phage therapy parameters

B = 1e6; % CFU, total bacteria
P = 0; % PFU, total phage
I = (2.7e6*lung_mass); % cells, initial number of neutrophils (from Roach et al., 2017)

% maximum number of neutrophils, we assume a large pool of neutrophils is available (10X Nlung)
nlung = 3.24e+06;
max_neutrophil_num = nlung*10; 

% bacterial inoculum and phage dose distribution
b_dist = 1:15; % uniformly distribute bacteria among 15 nodes
p_dist = 1:15; % uniformly distribute phage among 15 nodes

simu_time = 100; % simulation time in hours

% Simulate the metapopulation model
[time, res, p] = simu_metapop(B, P, I, max_neutrophil_num, b_dist, p_dist, simu_time);


%% Plotting metapopulation dynamics
figures_metapop(time,res,p)
