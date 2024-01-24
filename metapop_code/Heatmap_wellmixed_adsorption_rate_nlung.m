%% Analysis of the well-mixed model outcome given variations in immune levels and phage adsorption rate
% Inoculum: Phage-sensitive bacteria (B_S)
% Phage added two hours after infection

clc
clear;
close all; 

lung_mass = 0.135; % lung mass in grams

% Phage therapy parameters
B = 1e6; % CFU, total bacteria
P = 1e7; % PFU, total phage

% Number of lung neutrophils when host is immunocompetent
nlung = 3.24e+06;

% Percentage of neutrophil availability
perc = [0.01 0.1:0.05:1];

% Vector with phage adsorption rate values
adsorption_vec = linspace(-5, -9, 21);
adsorption_vec = 10.^adsorption_vec;

% total bacterial density
Tot_density = zeros(numel(adsorption_vec),  numel(perc));

simu_time = 250; % simulation time in hours

colors = brewermap(length(perc)+1, 'Paired');
colors = colors([1:length(perc)-1 length(perc)+1] ,:);

for ads = 1:numel(adsorption_vec)
    
    phage_ads = adsorption_vec(ads);
    
    for lvl = 1:numel(perc)   

        % Vary the number of neutrophil and initial immune density
        max_neutrophil_num = nlung*perc(lvl);
        if perc(lvl) >= 1
                 I = 2.7e6*lung_mass;
        else
                 I = (max_neutrophil_num/8.9); % initial amount of neutrophils is 8.9 times smaller than max_neutrophil_num
        end

        % Simulate the well-mixed model
        [time, res, p] = simu_metapop_singleNode_phageAdsorp(B, P, I, max_neutrophil_num, phage_ads, simu_time);
        
        BS_final = res(end, 1); % density
        BR_final = res(end, 16); % density
        num_bs = BS_final*sum(p.branch_volume.*p.nodes_pergen); % numbers
        num_br = BR_final*sum(p.branch_volume.*p.nodes_pergen); % numbers
        extinct_bs = num_bs < 1;
        extinct_br = num_br < 1;
        BS_final(extinct_bs) = 0;
        BR_final(extinct_br) = 0;

        Btot = BS_final + BR_final;
        Tot_density(ads, lvl) = Btot;

    end

end

save('../data/tot_bact_density_wellmixed_phageAdsorp.mat', 'Tot_density');

%%

load '../data/tot_bact_density_wellmixed_phageAdsorp.mat'
% Vector with phage adsorption rate values
adsorption_vec = linspace(-5, -9, 21);
adsorption_vec = 10.^adsorption_vec;
adsorption_vec = adsorption_vec(6:end);

Tot_density_copy = Tot_density(6:end,:);
Tot_density_copy(Tot_density_copy == 0) = 1;
Tot_density_copy = log10(Tot_density_copy);

% % heatmap
figure(1)
heat = imagesc(Tot_density_copy, [0 10]);
h = colorbar;
caxis([0 10])
cmap = colormap(parula(1e3));
cmap = [1,1,1; cmap];
colormap(cmap)
ylabel(h, 'log$_{10}$ Total Bacteria (CFU/ml)', 'interpreter', 'latex', 'fontsize', 17)


immune_ticks =  string(perc.*100) + '%';
xlabel('$\%$ of neutrophil availability in the lungs', 'interpreter', 'latex')
tot_xticks = length(immune_ticks);
set(gca,'XTick', [2:2:tot_xticks], 'XTickLabel', immune_ticks([2:2:tot_xticks]))
xtickangle(45)
adsorption_ticks = '10^{' + string(log10(adsorption_vec)) + '}';
ylabel('Phage adsorption rate ($(ml/PFU)^\sigma\, h^{-1}$) ', 'interpreter', 'latex')
set(gca,'YTick', 1:5:length(adsorption_vec), 'yticklabel', adsorption_ticks(1:5:length(adsorption_ticks)));

title({'Effects of varying $\widetilde{\phi}$ and innate immune levels'; 'on the well-mixed model outcome'}, 'FontSize', 18, 'interpreter', 'latex')
set(gca, 'fontsize', 17, 'linewidth', 1.5, 'TickDir','out')

% filename = '../figures/FigS_heatmap_well-mixed_phivsnlung.eps';
% exportgraphics(gcf, filename);
