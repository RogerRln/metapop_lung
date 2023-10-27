%% Read and analyze the anatomical data of the mouse bronchial tree (Counter et al., 2013)
% Extract the diameter and length of a branch for a specific generation
% Calculate L/D, D/DP, L/LP and the area of each branch at each generation
% using Counter et al. data

% Total bacteria in mice  lungs
W_lung = 0.173; % mice lungs weight
% generations = 13+1; % + 1 accounts for generation 0, the trachea
% w_generations = W_lung/generations; % weight per generation
% num_comp = 2.^(0:generations-1); % start at generation 0 (trachea), number of compartments per generation
% w_compartment = w_generations./num_comp; % weight per compartment at specific generation



carrying_cap = 10^10; % cfu/g, cell host &Microbe 2017 paper
total_bact = W_lung*carrying_cap; % total number of bacteria in mice lung


% Divide the original Kc by 14 generations, kc_uniformgen = carrying_cap/14;
% then at each generation the kc increases by 1/14 compare to previous
% generation
% at the end you have to compensate for not having a complete kc = 1e10 (which is needed
% for sum( weights_generations .* 1e10 ) = 1.73e9) per generation.
%, so at the last generation you add up all the missing proportions. 
% 
% kc_uniformgen = carrying_cap/100;% carrying capacity uniformly distributed in N generations
% 
% test = (kc_uniformgen + kc_uniformgen.*(0:0.1:1.2))
% sum(test)
%kc_pergen = kc_uniformgen.*(1:14);
%kc_pergen(end) = kc_pergen(end) + sum(kc_uniformgen.*(13:-1:1));

CT_Counter = readtable('./Counter_2012_airways.csv');
CT_Counter = CT_Counter(2:end,1:5);
CT_Counter = table2array(CT_Counter);
[r,c] = size(CT_Counter);
anatomy_airways = zeros(r,c);
for i = 1:r
    for j = 1:c
        strr = CT_Counter(i,j);
        if strr{1}
            anatomy_airways(i,j) = str2num(strr{1});
        else
            anatomy_airways(i,j) = 0;
        end
    end
end

uncertainty = @(x, y, SDx, SDy) abs(x./y).*sqrt((SDx./x).^2 + (SDy./y).^2);

generations = anatomy_airways(:,1);
diameter=  anatomy_airways(:,2);
diameter_std = anatomy_airways(:,3);
length = anatomy_airways(:,4);
length_std = anatomy_airways(:,5);

% Plot the anatomy changes of the airways through generations
% Diameter vs Generation
subplot(3,2,1);
errorbar(generations, diameter, diameter_std, 'o', 'color', 'k', 'linewidth', 1)
ylim([0 2])
xlim([0 16])
ylabel('Diameter (mm)')
xlabel('Generation')
set(gca, 'fontweight', 'bold', 'fontsize', 13)

% D/Dp vs Generation
subplot(3,2,3);
DDp_ratio =  diameter(2:end)./diameter(1:end-1);
DDpratio_std  = uncertainty(diameter(2:end), diameter(1:end-1), diameter_std(2:end), diameter_std(1:end-1));
errorbar(generations(2:end), DDp_ratio, DDpratio_std , 'o', 'color', 'k', 'linewidth', 1)
xlim([0 16])
ylim([0 2])
ylabel('D/D_P')
xlabel('Generation')
set(gca, 'fontweight', 'bold', 'fontsize', 13)


% L/D vs Generation
subplot(3,2,5);
LD_ratio = length./diameter;
LDratio_std  = uncertainty(length, diameter, length_std, diameter_std);
errorbar(generations, LD_ratio, LDratio_std , 'o', 'color', 'k', 'linewidth', 1)
xlim([0 16])
ylabel('L/D')
xlabel('Generation')
set(gca, 'fontweight', 'bold', 'fontsize', 13)

% L vs Generation
subplot(3,2,2);
errorbar(generations, length, length_std , 'o', 'color', 'k', 'linewidth', 1)
xlim([0 16])
ylabel('L (mm)')
xlabel('Generation')
set(gca, 'fontweight', 'bold', 'fontsize', 13)

% L/LP vs Generation
subplot(3,2,4);
LLP_ratio = length(2:end)./length(1:end-1);
LLPratio_std  = uncertainty(length(2:end), length(1:end-1), length_std(2:end), length_std(1:end-1));
errorbar(generations(2:end), LLP_ratio, LLPratio_std , 'o', 'color', 'k', 'linewidth', 1)
xlim([0 16])
ylabel('L/L_P')
xlabel('Generation')
set(gca, 'fontweight', 'bold', 'fontsize', 13)


% Area vs Generation
uncertainty_mul = @(x, y, SDx, SDy) abs(x.*y).*sqrt((SDx./x).^2 + (SDy./y).^2);
subplot(3,2,6);
Area = pi.*diameter.*length;
Area_std  = uncertainty_mul(pi.*diameter, length, pi*diameter_std, length_std);
errorbar(generations, Area, Area_std , 'o', 'color', 'k', 'linewidth', 1)
xlim([0 16])
ylabel('Area (mm^2)')
xlabel('Generation')
set(gca, 'fontweight', 'bold', 'fontsize', 13)
set(gcf,'position',[1000         692         800         646])

mucus_thickness = 10; %um
airways = (2.^(0:16))'; % number of airways per generation
total_mucusvolume = sum((Area*mucus_thickness*1e-3).*airways); % mm^3
%new_kc = total_bact/total_mucusvolume;
new_kc = total_bact/(total_mucusvolume*1e-3); %mm^3 to cm^3

% Calculate new carrying capacity based in Total surface area of the lung
% W_lung*carrying_cap = lung_surface*new_carrying_cap
lung_surface = sum(Area.*airways); % mm^2
%new_kc = total_bact/lung_surface; % cfu/mm^2

%% Scaling laws of diameter, length, area and volume 
% Diameter vs Generation
subplot(4,2,1);
errorbar(generations, diameter, diameter_std, 'o', 'color', 'k', 'linewidth', 1)
[params,S] = polyfit(generations, log(diameter), 1);
slope = params(1);
intercept = params(2);
exp_func = @(data) exp(slope.*data + intercept);
hold on
plot(generations, exp_func(generations), '--', 'color', 'r', 'linewidth', 1)
text( 8, 0.9, ['m = ' num2str(abs(params(1)))])
hold off
set(gca,'Yscale', 'log')
ylim([0 2])
xlim([0 16])
ylabel('Log Diameter (mm)')
xlabel('Generation')
set(gca, 'fontweight', 'bold', 'fontsize', 13)


% L vs Generation
subplot(4,2,3);
errorbar(generations, length, length_std , 'o', 'color', 'k', 'linewidth', 1)
[params,S] = polyfit(generations, log(length), 1);
slope = params(1);
intercept = params(2);
exp_func = @(data) exp(slope.*data + intercept);
hold on
plot(generations, exp_func(generations), '--', 'color', 'r', 'linewidth', 1)
text( 8, 7, ['m = ' num2str(abs(params(1)))])
hold off
set(gca,'Yscale', 'log')
xlim([0 16])
ylabel('Log L (mm)')
xlabel('Generation')
set(gca, 'fontweight', 'bold', 'fontsize', 13)

% Area vs Generation
uncertainty_mul = @(x, y, SDx, SDy) abs(x.*y).*sqrt((SDx./x).^2 + (SDy./y).^2);
subplot(4,2,5);
Area = pi.*diameter.*length;
Area_std  = uncertainty_mul(pi.*diameter, length, pi*diameter_std, length_std);
errorbar(generations, Area, Area_std , 'o', 'color', 'k', 'linewidth', 1)
[params,S] = polyfit(generations, log(Area), 1);
slope = params(1);
intercept = params(2);
exp_func = @(data) exp(slope.*data + intercept);
hold on
plot(generations, exp_func(generations), '--', 'color', 'r', 'linewidth', 1)
text( 8, 20, ['m = ' num2str(abs(params(1)))])
hold off
set(gca,'Yscale', 'log')
xlim([0 16])
ylabel('Log Area (mm^2)')
xlabel('Generation')
set(gca, 'fontweight', 'bold', 'fontsize', 13)

% Plot the anatomy changes of the airways through generations
% Diameter vs Generation
subplot(4,2,2);
errorbar(generations, diameter, diameter_std, 'o', 'color', 'k', 'linewidth', 1)
f = fit(generations(2:end), diameter(2:end), 'power1'); % fit a power law of the form a*x^b
base = f.a;
exponent = f.b;
fit_data= base.*generations(2:end).^exponent;
results.base = base;
hold on
plot(generations(2:end), fit_data, 'r--')
text( 8, 1.1, ['exp = ' num2str(abs(exponent))])
hold off
set(gca,'Yscale', 'log')
set(gca,'Xscale', 'log')
ylim([0 2])
xlim([0 16])
ylabel('Log Diameter (mm)')
xlabel('Log Generation')
set(gca, 'fontweight', 'bold', 'fontsize', 13)


% L vs Generation
subplot(4,2,4);
errorbar(generations, length, length_std , 'o', 'color', 'k', 'linewidth', 1)
f = fit(generations(2:end), length(2:end), 'power1'); % fit a power law of the form a*x^b
base = f.a;
exponent = f.b;
fit_data= base.*generations(2:end).^exponent;
results.base = base;
hold on
plot(generations(2:end), fit_data, 'r--')
text( 8, 6, ['exp = ' num2str(abs(exponent))])
hold off
set(gca,'Yscale', 'log')
set(gca,'Xscale', 'log')
xlim([0 16])
ylabel('Log L (mm)')
xlabel('Log Generation')
set(gca, 'fontweight', 'bold', 'fontsize', 13)

% Area vs Generation
uncertainty_mul = @(x, y, SDx, SDy) abs(x.*y).*sqrt((SDx./x).^2 + (SDy./y).^2);
subplot(4,2,6);
Area = pi.*diameter.*length;
Area_std  = uncertainty_mul(pi.*diameter, length, pi*diameter_std, length_std);
errorbar(generations, Area, Area_std , 'o', 'color', 'k', 'linewidth', 1)
f = fit(generations(2:end), Area(2:end), 'power1'); % fit a power law of the form a*x^b
base = f.a;
exponent = f.b;
fit_data= base.*generations(2:end).^exponent;
hold on
plot(generations(2:end), fit_data, 'r--')
text( 8, 12, ['exp = ' num2str(abs(exponent))])
hold off
set(gca,'Yscale', 'log')
set(gca,'Xscale', 'log')
xlim([0 16])
ylabel('Log Area (mm^2)')
xlabel('Log Generation')
set(gca, 'fontweight', 'bold', 'fontsize', 13)


% Volume vs Generation
uncertainty_mul = @(x, y, SDx, SDy) abs(x.*y).*sqrt((SDx./x).^2 + (SDy./y).^2);

radius = diameter./2;
radius_square = radius.^2;
radius_std = diameter_std./2;
n = 2;
uncertainty_radiusquare = @(r) abs(n).*(radius_std./abs(r)).*abs(radius_square);

Volume = pi.*radius_square.*length;
Volume_std  = uncertainty_mul(pi.*radius_square, length, pi*uncertainty_radiusquare(radius), length_std);

subplot(4,2,7);
errorbar(generations, Volume, Volume_std , 'o', 'color', 'k', 'linewidth', 1)
[params,S] = polyfit(generations, log(Volume), 1);
slope = params(1);
intercept = params(2);
exp_func = @(data) exp(slope.*data + intercept);
hold on
plot(generations, exp_func(generations), '--', 'color', 'r', 'linewidth', 1)
text( 8, 3, ['m = ' num2str(abs(params(1)))])
hold off
set(gca,'Yscale', 'log')
xlim([0 16])
ylabel('Log Vol (mm^3)')
xlabel('Generation')
set(gca, 'fontweight', 'bold', 'fontsize', 13)



% Volume vs Generation
uncertainty_mul = @(x, y, SDx, SDy) abs(x.*y).*sqrt((SDx./x).^2 + (SDy./y).^2);

radius = diameter./2;
radius_square = radius.^2;
radius_std = diameter_std./2;
n = 2;
uncertainty_radiusquare = @(r) abs(n).*(radius_std./abs(r)).*abs(radius_square);

Volume = pi.*radius_square.*length;
Volume_std  = uncertainty_mul(pi.*radius_square, length, pi*uncertainty_radiusquare(radius), length_std);

subplot(4,2,8);

errorbar(generations, Volume, Volume_std , 'o', 'color', 'k', 'linewidth', 1)
f = fit(generations(2:end), Volume(2:end), 'power1'); % fit a power law of the form a*x^b
base = f.a;
exponent = f.b;
fit_data= base.*generations(2:end).^exponent;
hold on
plot(generations(2:end), fit_data, 'r--')
text( 7, 4, ['exp = ' num2str(abs(exponent))])
hold off
set(gca,'Yscale', 'log')
set(gca,'Xscale', 'log')
xlim([0 16])
ylabel('Log vol (mm^3)')
xlabel('Log Generation')
set(gca, 'fontweight', 'bold', 'fontsize', 13)



set(gcf,'position',[1000         692         800         646])


%% Anatomical measurements per generation (Lenght, Diameter, Area, Volume)



% Plot the anatomy changes of the airways through generations
% Diameter vs Generation
subplot(2,2,1);
errorbar(generations, diameter, diameter_std, 'o', 'color', 'k', 'linewidth', 1)
ylim([0 2])
xlim([0 16])
ylabel('Diameter (mm)')
xlabel('Generation')
set(gca, 'fontweight', 'bold', 'fontsize', 13)

% Length vs Generation
subplot(2,2,2);
errorbar(generations, length, length_std , 'o', 'color', 'k', 'linewidth', 1)
xlim([0 16])
ylabel('L (mm)')
xlabel('Generation')
set(gca, 'fontweight', 'bold', 'fontsize', 13)



% Area vs Generation
uncertainty_mul = @(x, y, SDx, SDy) abs(x.*y).*sqrt((SDx./x).^2 + (SDy./y).^2);
subplot(2,2,3);
Area = pi.*diameter.*length;
Area_std  = uncertainty_mul(pi.*diameter, length, pi*diameter_std, length_std);
errorbar(generations, Area, Area_std , 'o', 'color', 'k', 'linewidth', 1)
xlim([0 16])
ylabel('Area (mm^2)')
xlabel('Generation')
set(gca, 'fontweight', 'bold', 'fontsize', 13)


% Volume vs Generation, 
uncertainty_mul = @(x, y, SDx, SDy) abs(x.*y).*sqrt((SDx./x).^2 + (SDy./y).^2);
radius = diameter./2;
radius_square = radius.^2;
radius_std = diameter_std./2;
n = 2;
uncertainty_radiusquare = @(r) abs(n).*(radius_std./abs(r)).*abs(radius_square);
Volume = pi.*radius_square.*length;
Volume_std  = uncertainty_mul(pi.*radius_square, length, pi*uncertainty_radiusquare(radius), length_std);

subplot(2,2,4);
errorbar(generations, Volume, Volume_std , 'o', 'color', 'k', 'linewidth', 1)
xlim([0 16])
ylabel('Volume (mm^3)')
xlabel('Generation')
set(gca, 'fontweight', 'bold', 'fontsize', 13)
set(gcf,'position',[1000         692         800         646])

%% Scaling laws volume


% Volume vs Generation, 
uncertainty_mul = @(x, y, SDx, SDy) abs(x.*y).*sqrt((SDx./x).^2 + (SDy./y).^2);

radius = diameter./2;
radius_square = radius.^2;
radius_std = diameter_std./2;
n = 2;
uncertainty_radiusquare = @(r) abs(n).*(radius_std./abs(r)).*abs(radius_square);

Volume = pi.*radius_square.*length;
Volume_std  = uncertainty_mul(pi.*radius_square, length, pi*uncertainty_radiusquare(radius), length_std);

subplot(3,1,1);
errorbar(generations, Volume, Volume_std , 'o', 'color', 'k', 'linewidth', 1)
xlim([0 16])
ylabel('Volume (mm^3)')
xlabel('Generation')
set(gca, 'fontweight', 'bold', 'fontsize', 13)

% LOG Volume vs Generation
uncertainty_mul = @(x, y, SDx, SDy) abs(x.*y).*sqrt((SDx./x).^2 + (SDy./y).^2);

radius = diameter./2;
radius_square = radius.^2;
radius_std = diameter_std./2;
n = 2;
uncertainty_radiusquare = @(r) abs(n).*(radius_std./abs(r)).*abs(radius_square);

Volume = pi.*radius_square.*length;
Volume_std  = uncertainty_mul(pi.*radius_square, length, pi*uncertainty_radiusquare(radius), length_std);

subplot(3,1,2);
errorbar(generations, Volume, Volume_std , 'o', 'color', 'k', 'linewidth', 1)
[params,S] = polyfit(generations, log(Volume), 1);
slope = params(1);
intercept = params(2);
exp_func = @(data) exp(slope.*data + intercept);
hold on
plot(generations, exp_func(generations), '--', 'color', 'r', 'linewidth', 1)
text( 8, 3, ['m = ' num2str(abs(params(1)))])
set(gca,'Yscale', 'log')
xlim([0 16])
ylabel('Log vol (mm^3)')
xlabel('Generation')
set(gca, 'fontweight', 'bold', 'fontsize', 13)



% LOG Volume vs LOG Generation
uncertainty_mul = @(x, y, SDx, SDy) abs(x.*y).*sqrt((SDx./x).^2 + (SDy./y).^2);

radius = diameter./2;
radius_square = radius.^2;
radius_std = diameter_std./2;
n = 2;
uncertainty_radiusquare = @(r) abs(n).*(radius_std./abs(r)).*abs(radius_square);

Volume = pi.*radius_square.*length;
Volume_std  = uncertainty_mul(pi.*radius_square, length, pi*uncertainty_radiusquare(radius), length_std);

subplot(3,1,3);

errorbar(generations, Volume, Volume_std , 'o', 'color', 'k', 'linewidth', 1)
f = fit(generations(2:end), Volume(2:end), 'power1'); % fit a power law of the form a*x^b
base = f.a;
exponent = f.b;
fit_data= base.*generations(2:end).^exponent;
hold on
plot(generations(2:end), fit_data, 'r--')
text( 7, 4, ['exp = ' num2str(abs(exponent))])
hold off
set(gca,'Yscale', 'log')
set(gca,'Xscale', 'log')
xlim([0 16])
ylabel('Log vol (mm^3)')
xlabel('Log Generation')
set(gca, 'fontweight', 'bold', 'fontsize', 13)
set(gcf,'position',[1000         692         800         646])