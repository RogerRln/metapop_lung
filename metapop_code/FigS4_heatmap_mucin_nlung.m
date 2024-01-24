%% Code to perform the robustness analysis of the model by varying mucin and immune levels
% Inoculum: Phage-sensitive bacteria (B_S)
% Phage added two hours after infection

% If you want to run everything from scratch uncomment lines 40 to 69
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

% Total number of lung neutrophils
nlung = 3.24e+06;

% mucin concentration (%)
mucin_conc = linspace(4,0,17);

% Load phage and bacteria distribution parameters
load inocuDist_params.mat

simu_time = 250; % simulation time in hours

% Matrix  with prob. of therapeutic success due to variations in initial
%conditions and given a specific phage adsorption rate and immune level
therapeutic_success = zeros(length(mucin_conc), length(perc));

% Matrix with median time to bacterial extinction 
Time_Extinction = zeros(length(mucin_conc), length(perc));

% tic
% for mu = 1:numel(mucin_conc)   
%     
%     % variations in the phage adsorption rate
%     mucin = mucin_conc(mu);
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
%         [Total_CFU_params, Time_extinct_vec, p] = robust_mucin_immune(B, P, I, max_neutrophil_num, simu_time, mucin, inocuDist_params);
%         
%         % calculate probability of therapeutic success
%         prob_succ = sum(Total_CFU_params < 1)/numel(Total_CFU_params);
%         therapeutic_success(mu, lvl) = prob_succ;
%         % median time to bacterial extinction
%         Time_Extinction(mu, lvl) = median(Time_extinct_vec);   
%     end
% end
% toc

% save('Prob_success_mucin_vs_nlung_test.mat', 'therapeutic_success');
% save('Median_TimeExtinction_mucin_vs_nlung_test.mat', 'Time_Extinction');

%% Plotting probability of therapeutic success 

load '../data/Prob_success_mucin_vs_nlung_test.mat'
load '../data/Median_TimeExtinction_mucin_vs_nlung_test.mat'

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
xlabel('$\%$ of neutrophil availability in the lungs', 'interpreter', 'latex')
tot_xticks = length(immune_ticks);
set(gca,'XTick', [2:2:tot_xticks], 'XTickLabel', immune_ticks([2:2:tot_xticks]))
xtickangle(45)

ylabel('Mucin concentration (\%)', 'interpreter', 'latex')
set(gca,'YTick', 1:4:length(mucin_conc), 'yticklabel', mucin_conc(1:4:length(mucin_conc)));


title({'Effects of varying mucin and innate immune levels'; 'on the spatial model outcome'}, 'FontSize', 18, 'interpreter', 'latex')
set(gca, 'fontsize', 17, 'linewidth', 1.5, 'TickDir','out')


% filename = '/Users/rrodriguez77/Dropbox (GaTech)/Phage-Immune_host-pathogen project/draft_metapopulation/version2_Sep132022/figures/variation_initial_conditions/PhageAdsorp_vs_nlung/prob_success.eps';
% exportgraphics(gcf, filename);

% Time to extinction heatmap
figure(2)
time_ext_copy = Time_Extinction;
heat = imagesc(time_ext_copy);
h = colorbar;
caxis([min(min(Time_Extinction))-5 max(max(time_ext_copy))])
%caxis([min(min(Time_Extinction)) max(max(Time_Extinction))])
cmap = colormap(parula(1e3));
colormap(cmap)
ylabel(h, 'Time (h)', 'interpreter', 'latex', 'fontsize', 14)

xlabel('$\%$ of neutrophil availability in the lungs', 'interpreter', 'latex')
tot_xticks = length(immune_ticks);
set(gca,'XTick', [2:2:tot_xticks], 'XTickLabel', immune_ticks([2:2:tot_xticks]))
xtickangle(45)

ylabel('Mucin concentration (%)', 'interpreter', 'latex', 'fontsize', 15)
set(gca,'YTick', 1:4:length(mucin_conc), 'yticklabel', mucin_conc(1:4:length(mucin_conc)));

title('Median of time to extinction', 'FontSize', 14, 'interpreter', 'latex')
set(gca, 'fontsize', 15, 'linewidth', 1.5)
hold on
plot([9 9], [0 numel(mucin_conc)+1], 'linewidth', 1.5, 'color', 'black');
hold off

% filename = '/Users/rrodriguez77/Dropbox (GaTech)/Phage-Immune_host-pathogen project/draft_metapopulation/version2_Sep132022/figures/variation_initial_conditions/median_time_extinction.eps';
% exportgraphics(gcf, filename);