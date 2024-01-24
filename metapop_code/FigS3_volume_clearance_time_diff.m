%% Code to calculate bacterial extinction times due to phage therapy of P.a. infection in an immunocompetent host

clc
clear;
close all; 

%% We first simulate the phage therapy of a P.a. infection in an immunocompetent host
% Inoculum: Phage-sensitive bacteria (B_S)
% Phage added two hours after infection

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

%% Calculate bacterial extinction times from simulations and from phage lysis and immune killing rates

eli_time = time_clearance_thresh(res, time, p, 1./p.branch_volume); % elimination times from simulations
log_volratio = log(p.branch_volume(1:end-1)./p.branch_volume(end)); % natural log of volume ratio V_i/V_bottom
time_diff = eli_time(1:end-1) - eli_time(end); % elimination time difference between node i and terminal node

bs = res(1:end, 1:15); % phage-sensitive bacteria density
ph = res(1:end, (2*p.NP+1):3*p.NP); % phage density


% calculate immune killing rate (maximum immune killing at the end of the simulation)
immune_killing = p.ep*p.Ki;


%% Plot infection clearance time differences between node i nad terminal node

figure(1);
plot(log_volratio, time_diff , 'ko', 'LineWidth', 2.5)
parms = polyfit(log_volratio, time_diff, 1);
hold on
plot([6;log_volratio;0], (1/(immune_killing - p.rs)).*[6;log_volratio;0], '--', 'LineWidth', 2.5, 'color', [0.4940 0.1840 0.5560])
hold off
ylim([0 5])
xlabel('ln($\frac{V_i}{V_{bottom}}$)', 'Interpreter', 'latex')
ylabel('$T_i - T_{bottom}$ (h)', 'Interpreter', 'latex')
legend('Simulations', 'Immune killing', 'location', 'northwest')
legend box off
set(gca, 'fontsize', 18, 'linewidth', 1.5)
set(gcf, 'renderer', 'painters')
% file_name = '/Users/rrodriguez77/Dropbox (GaTech)/Phage-Immune_host-pathogen project/draft_metapopulation/version2_Sep132022/figures/model_volume_theory/volume_vs_clearance_node3_t41_densities.eps';
% exportgraphics(gcf,file_name);

%% Find the simulation time that the minimize the error between theoretical elimination time
% due to phage lysis + immune killing and the elimination time obtained directly from the simulations

% bs = res(1:end, 1:15);
% br = res(1:end, 16:30);
% ph = res(1:end, (2*p.NP+1):3*p.NP);
% 
% node_id = 5;
% time_id = 1;
% 
% for j = 1:numel(bs(:,1))
%     phage_killing = (p.phi*(ph(j, node_id))^p.g);
%     immune_killing = (p.ep*p.Ki)/(1 + (bs(j, node_id))/p.Kd);
%     both_killing = phage_killing + immune_killing;
%     theory_rate  = (1/both_killing).*log_volratio;
%     mse(j) = mean((time_diff-theory_rate).^2);
% end
% 
% [c, ind] = sort(mse);
% 
% 
% phage_killing = (p.phi*(ph(ind(time_id), node_id))^p.g);
% immune_killing = (p.ep*p.Ki)/(1+ (bs(ind(time_id), node_id))/p.Kd);
% both_killing = phage_killing + immune_killing;
% 
% 
% nid = 3;
% tid = 41;
%nid = 1;
%tid = 26; 
%nid = 3;
%tid = 35;