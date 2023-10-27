%% Code to simulate the well-mixed (single node) model given variations in immune levels
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
perc = [0.1:0.1:0.4 0.45 0.5:0.1:1];

simu_time = 115; % simulation time in hours

colors = brewermap(length(perc)+1, 'Paired');
colors = colors([1:length(perc)-1 length(perc)+1] ,:);

for lvl = 1:numel(perc)   
        
    % Vary the number of neutrophil and initial immune density
    max_neutrophil_num = nlung*perc(lvl);
    if perc(lvl) >= 1
             I = 2.7e6*lung_mass;
    else
             I = (max_neutrophil_num/8.9); % initial amount of neutrophils is 8.9 times smaller than max_neutrophil_num
    end

    % Simulate the well-mixed model
    [time, res, p] = simu_metapop_singleNode(B, P, I, max_neutrophil_num, simu_time);
    Btot = sum(res(:, [1 16]), 2); % B_S + B_R
    
    % Plot total bacterial density
    semilogy(time, Btot, 'linewidth', 2.5, 'color', colors(lvl,:))
    hold on
    ylim([1/p.lung_volume 1e10])
    xlim([0 simu_time])
    xlabel('Time (h)', 'interpreter', 'latex')
    ylabel('CFU/ml', 'interpreter', 'latex')
end

hold off
labels =  string(perc*100) + '\% of $N_{Lung}$';
legend(labels, 'location', 'best', 'interpreter', 'latex')
legend boxoff
set(gca, 'fontsize', 15, 'linewidth', 1.5)


