
% low mucin conc. [0-1%]

mucin = [0	0.2	0.6	1];
%mucin = [0	0.2/100	0.6/100	1/100];
visc_M = [0.69	0.80	0.96	0.99];
visc_W = 0.69;

y = log(visc_M-visc_W);
x = log(mucin);

[p,S] = polyfit(x(2:end),y(2:end),1);
a = exp(p(2)); % intercept
b = p(1); % slope

y1 = polyval(p, x(2:end));

% exp(y1) = (visc_M - visc_W)
% visc_M = exp(y1) + visc_W
% visc_M = a*mucin^b + visc_W

pred_visc_M = a.*mucin.^b + visc_W;

% high mucin conc. [1-4%]

mucin = [1 2 4];
visc_M = [0.99	1.83	3.39];
visc_W = 0.69;

y = log(visc_M-visc_W);
x = log(mucin);

[p,S] = polyfit(x,y,1);
a = exp(p(2)); % intercept
b = p(1); % slope

y1 = polyval(p, x);

% exp(y1) = (visc_M - visc_W)
% visc_M = exp(y1) + visc_W
% visc_M = a*mucin^b + visc_W

pred_visc_M = a.*mucin.^b + visc_W;

high_visc = @(mu) a.*mu.^b + visc_W;

30/(6*pi*0.68*high_visc(2.5))

30/(6*pi*0.68*high_visc(8))

%% plot mucin concentration vs Diffusion constant without datapoints

load struct_variable.mat
p = struct_variable;

mucin_conc = 0:0.1:8;
%mucin_conc = [0 0.2 0.6 1 2 4];

BD = zeros(size(mucin_conc));
PD = zeros(size(mucin_conc));
ND = zeros(size(mucin_conc));

for i = 1:length(mucin_conc)
    [bac_diff, phage_d, n_speed] = mucin_to_Diff(mucin_conc(i));
    BD(i) = bac_diff;
    PD(i) = phage_d;
    speed_neutro = n_speed*1e-4*60; % from um/min to cm/h
    n_diff = (speed_neutro^2)/(2*p.alpha);
    n_diff = (n_diff*1e8)/3600; % from cm^2/h to um^2/sec
    ND(i) = n_diff;
end

cols = [0 0.4470 0.7410; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560];
figure(7)
semilogy(mucin_conc, [BD; PD; ND], 'linewidth', 2.5)% I am playing to see if i recover neutrophil speed from diffusion
title('Mucin concentration vs Diffusion constant')
legend('Bacteria', 'Phage', 'Neutrophils' , 'Location', 'best')
legend box off
set(gca, 'fontsize', 15, 'Linewidth', 1.5, 'ColorOrder', cols)
xlabel('Mucin concentration (%)')
ylabel('Diffusion constant ($\mu m^2/s$)', 'interpreter', 'latex')

%% plot mucin concentration vs Diffusion constant WITH data points


mucin_conc = 0:0.1:8;
%mucin_conc = [0 0.2 0.6 1 2 4];

BD = zeros(size(mucin_conc));
PD = zeros(size(mucin_conc));
ND = zeros(size(mucin_conc));

for i = 1:length(mucin_conc)
    [bac_diff, phage_d, n_speed] = mucin_to_Diff(mucin_conc(i));
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
mucin_neutro = [1.5 2.5 6.5];
neutro_D = zeros(numel(mucin_neutro),1);
for i = 1:numel(mucin_neutro)
    [~, ~, n_speed] = mucin_to_Diff(mucin_neutro(i));
    speed_neutro = n_speed*1e-4*60; % from um/min to cm/h
    n_diff = (speed_neutro^2)/(2*p.alpha);
    n_diff = (n_diff*1e8)/3600; % from cm^2/h to um^2/sec
    neutro_D(i) = n_diff;
end


cols = [0 0.4470 0.7410; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560];
figure(1)
%semilogy(mucin_conc, [BD; PD; ND], 'linewidth', 2.5)
semilogy(mucin_conc, [BD; PD], 'linewidth', 2.5)
hold on
semilogy(mucin_bac, bac_D, 'ok', 'linewidth', 2, 'HandleVisibility', 'off', 'MarkerFaceColor', 'k')
semilogy(mucin_phage, phage_D, 'ok', 'linewidth', 2, 'MarkerFaceColor', 'k')
%semilogy(mucin_neutro, neutro_D, 'ok', 'linewidth', 2, 'MarkerFaceColor', 'k')
rectangle('Position',[4.1,1e-1,3.9,1e3],'FaceColor',[128/255 128/255 128/255 0.4],'LineWidth',1.5)
hold off
title('Mucin level vs Diffusion constant','interpreter', 'latex', 'fontsize', 19)
%legend('Bacteria', 'Phage', 'Neutrophils', 'Data', 'Location', 'best')
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
bact_speed = [17 1]; % units um/sec
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


figure(3)
mucin = [1.5 2.5 6.5];
neutro_speed = [5 3 2.3]; % units um/min

query_high = 2.5:0.5:6.5;
interp_speed_high = interp1(mucin, neutro_speed, query_high);
[par_neutro_high,S] = polyfit(query_high,interp_speed_high,1);

query_low = 1.5:0.5:2.5;
interp_speed_low = interp1(mucin, neutro_speed, query_low);
[par_neutro_low,S] = polyfit(query_low,interp_speed_low,1);

n_speed_low = @(mu) par_neutro_low(1).*mu + par_neutro_low(2);
n_speed_high = @(mu) par_neutro_high(1).*mu + par_neutro_high(2);


plot(mucin_conc, [n_speed_low(0:0.1:2.4) n_speed_high(2.5:0.1:8)], 'color', [0.4940 0.1840 0.5560], 'linewidth', 2.5)
hold on
plot(query_low, interp_speed_low, 'r--', 'LineWidth', 2)
plot(query_high, interp_speed_high, 'g--', 'linewidth', 2)
plot(mucin, neutro_speed, 'ok', 'linewidth', 2, 'MarkerFaceColor', 'k')
hold off
title('Mucin concentration vs neutrophil speed')
legend('Neutrophil speed', 'fitted line 1', 'fitted line 2', 'Data', 'Location', 'best')
legend box off
set(gca, 'fontsize', 15, 'Linewidth', 1.5, 'ColorOrder', cols)
xlabel('Mucin concentration (%)')
ylabel('Speed ($\mu m/min$)', 'interpreter', 'latex')
% filename = '/Users/rrodriguez77/Dropbox (GaTech)/Phage-Immune_host-pathogen project/draft_metapopulation/version2_Sep132022/figures/Mucin_vs_Diff/mucin_vs_NeutroSpeed.eps';
% exportgraphics(gcf, filename);