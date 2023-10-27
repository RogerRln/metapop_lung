%% Code to perform the robustness analysis of the model
% We vary innate immune levels and phage adsorption rate and calculate the
% probability of therapeutic success
% Inoculum: Phage-sensitive bacteria (B_S)
% Phage added two hours after infection

% If you want to run everything from scratch uncomment lines 45 to 73
% it will take ~13 hours to perform the robustness analysis

clc
clear;
close all; 


lung_mass = 0.135; % lung mass in grams

% Initial bacteria and phage dose
B = 1e6; % CFU, total bacteria
P = 1e7; % PFU, total phage

% Percentage of neutrophil availability
perc = [0.01 0.1:0.05:1];

% total number of lung neutrophils
nlung = 3.24e+06;

% Vector with phage adsorption rate values
adsorption_vec = linspace(-5, -9, 21);
adsorption_vec = 10.^adsorption_vec;

% Load phage and bacteria distribution parameters
load inocuDist_params.mat


simu_time = 250; % simulation time in hours

% Matrix  with prob. of therapeutic success due to variations in initial
%conditions and given a specific phage adsorption rate and immune level
therapeutic_success = zeros(length(adsorption_vec), length(perc));

% Matrix with median time to bacterial extinction 
Time_Extinction = zeros(length(adsorption_vec), length(perc));

% tic
% for ads = 1:numel(adsorption_vec)   
%     
%     % variations in the phage adsorption rate
%     phage_ads = adsorption_vec(ads);
%     
%     for lvl = 1:numel(perc)
%         
%         % Varying the number of neutrophils and initial immune density
%         max_neutrophil_num = nlung*perc(lvl);
%         if perc(lvl) >= 1
%                  I = 2.7e6*lung_mass;
%         else
%                  I = (max_neutrophil_num/8.9); % initial amount of neutrophils is 8.9 times smaller than max_neutrophil_num
%         end
%  
%         % Simulate the model
%         [Total_CFU_params, Time_extinct_vec, p] = robust_adsorption_immune(B, P, I, max_neutrophil_num, simu_time, phage_ads, inocuDist_params);
%         
%         % calculate probability of therapeutic success
%         prob_succ = sum(Total_CFU_params < 1)/numel(Total_CFU_params);
%         therapeutic_success(ads, lvl) = prob_succ;
%         % median time to bacterial extinction
%         Time_Extinction(ads, lvl) = median(Time_extinct_vec);   
%     end
% end
% toc
% save('Prob_success_PhageAdsorp_vs_nlung_test.mat', 'therapeutic_success');
% save('Median_TimeExtinction_PhageAdsorp_vs_nlung_test.mat', 'Time_Extinction');

%% Plotting probability of therapeutic success 

load '../data/Prob_success_PhageAdsorp_vs_nlung_test.mat'

% heatmap
figure(1)
heat = imagesc(therapeutic_success, [0 1]);
h = colorbar;
caxis([0 1])
cmap = colormap(parula(1e3));
cmap = [0,0,0; cmap];
colormap(cmap)
ylabel(h, 'Probability of therapeutic success', 'interpreter', 'latex', 'fontsize', 17)

immune_ticks =  string(perc.*100) + '%';
xlabel('$N_{lung}$', 'interpreter', 'latex')
tot_xticks = length(immune_ticks);
set(gca,'XTick', [2:2:tot_xticks], 'XTickLabel', immune_ticks([2:2:tot_xticks]))
xtickangle(45)
adsorption_ticks = '10^{' + string(log10(adsorption_vec)) + '}';
ylabel('Phage adsorption rate ($(ml/PFU)^\sigma\, h^{-1}$) ', 'interpreter', 'latex')
set(gca,'YTick', 1:5:length(adsorption_vec), 'yticklabel', adsorption_ticks(1:5:length(adsorption_ticks)));

title('Probability of success in clearing the infection', 'FontSize', 18, 'interpreter', 'latex')
set(gca, 'fontsize', 17, 'linewidth', 1.5, 'TickDir','out')

% filename = '/Users/rrodriguez77/Dropbox (GaTech)/Phage-Immune_host-pathogen project/draft_metapopulation/version2_Sep132022/figures/variation_initial_conditions/PhageAdsorp_vs_nlung/prob_success.eps';
% exportgraphics(gcf, filename);