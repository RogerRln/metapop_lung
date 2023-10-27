%% Code to simulate the phage therapy of a P.a. infection in an immunocompetent host
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

% Simulate the metapopulation model
[time, res, p] = simu_metapop(B, P, I, max_neutrophil_num, b_dist, p_dist, simu_time);


%% Plotting metapopulation dynamics
figures_metapop(time,res,p)
