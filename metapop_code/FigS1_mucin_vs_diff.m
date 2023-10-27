%% Code to calculate and plot the diffusion constants of phage and bacteria
% as a function of mucin level

clc
clear;
close all;


%% plot mucin concentration vs Diffusion constant with observed data points

load struct_variable.mat
p = struct_variable;

mucin_conc = 0:0.1:8; % mucin concentration 0-8%

BD = zeros(size(mucin_conc));
PD = zeros(size(mucin_conc));
ND = zeros(size(mucin_conc));

for i = 1:length(mucin_conc)
    [bac_diff, phage_d, n_speed] = mucin_to_Diff(mucin_conc(i)); % mucin level to diffusion
    BD(i) = bac_diff;
    PD(i) = phage_d;
    speed_neutro = n_speed*1e-4*60; % from um/min to cm/h
    n_diff = (speed_neutro^2)/(2*p.alpha);
    n_diff = (n_diff*1e8)/3600; % from cm^2/h to um^2/sec
    ND(i) = n_diff;
end

% Obtain observed diffusion values for bacteria 
mucin_bac = [2.5 8];
bac_D = zeros(numel(mucin_bac),1);
for i = 1:numel(mucin_bac)
    [bac_diff, ~, ~] = mucin_to_Diff(mucin_bac(i));
    bac_D(i) = bac_diff;
end

% Obtain observed diffusion values for bacteria 
mucin_phage = [0 0.2 0.6 1 2 4];
phage_D = zeros(numel(mucin_phage),1);
for i = 1:numel(mucin_phage)
    [~, phage_d, ~] = mucin_to_Diff(mucin_phage(i));
    phage_D(i) = phage_d;
end

% Obtain observed diffusion values for neutrophils 
% mucin_neutro = [1.5 2.5 6.5];
% neutro_D = zeros(numel(mucin_neutro),1);
% for i = 1:numel(mucin_neutro)
%     [~, ~, n_speed] = mucin_to_Diff(mucin_neutro(i));
%     speed_neutro = n_speed*1e-4*60; % from um/min to cm/h
%     n_diff = (speed_neutro^2)/(2*p.alpha);
%     n_diff = (n_diff*1e8)/3600; % from cm^2/h to um^2/sec
%     neutro_D(i) = n_diff;
% end


cols = [0 0.4470 0.7410; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560];

figure(1)
semilogy(mucin_conc, [BD; PD], 'linewidth', 2.5)
hold on
semilogy(mucin_bac, bac_D, 'ok', 'linewidth', 2, 'HandleVisibility', 'off', 'MarkerFaceColor', 'k')
semilogy(mucin_phage, phage_D, 'ok', 'linewidth', 2, 'MarkerFaceColor', 'k')
rectangle('Position',[4.1,1e-1,3.9,1e3],'FaceColor',[128/255 128/255 128/255 0.4],'LineWidth',1.5)
hold off
title('Mucin level vs Diffusion constant','interpreter', 'latex', 'fontsize', 19)
legend('Bacteria', 'Phage', 'Data', 'Location', 'best')
legend box off
set(gca, 'fontsize', 17, 'Linewidth', 1.5, 'ColorOrder', cols)
xlabel('Mucin concentration (\%)', 'interpreter', 'latex')
ylabel('Diffusion constant ($\mu m^2/s$)', 'interpreter', 'latex')
set(gcf, 'renderer', 'painters')
% filename = '/Users/rrodriguez77/Dropbox (GaTech)/Phage-Immune_host-pathogen project/draft_metapopulation/version2_Sep132022/figures/Mucin_vs_Diff/mucin_vs_diff_wData.png';
% exportgraphics(gcf, filename, 'Resolution', 300);


figure(2)
mucin = [2.5 8];
bact_speed = [17 1]; % units um/sec, from Matsui H. et al.,(2006). PNAS
query = 2.5:0.1:8;
interp_speed = interp1(mucin, bact_speed, query);
[par_bac,S] = polyfit(query, interp_speed, 1);
speed = @(mu) par_bac(1).*mu + par_bac(2);

plot(mucin_conc, speed(mucin_conc), 'color', [0 0.4470 0.7410], 'linewidth', 2.5)
hold on
plot(query, interp_speed, 'r--', 'LineWidth', 2)
plot(mucin, bact_speed, 'ok', 'linewidth', 2, 'MarkerFaceColor', 'k')
rectangle('Position',[4.1,0,3.9,25],'FaceColor',[128/255 128/255 128/255 0.4],'LineWidth',1.5)
hold off
title('Mucin level vs Bacterial speed', 'interpreter', 'latex', 'fontsize', 19)
legend('Bacterial speed', 'fitted line', 'Data', 'Location', 'best')
legend box off
set(gca, 'fontsize', 17, 'Linewidth', 1.5, 'ColorOrder', cols)
xlabel('Mucin concentration (\%)', 'interpreter', 'latex')
ylabel('Speed ($\mu m/s$)', 'interpreter', 'latex')
set(gcf, 'renderer', 'painters')
% filename = '/Users/rrodriguez77/Dropbox (GaTech)/Phage-Immune_host-pathogen project/draft_metapopulation/version2_Sep132022/figures/Mucin_vs_Diff/mucin_vs_BactSpeed.png';
% exportgraphics(gcf, filename, 'Resolution',300);

