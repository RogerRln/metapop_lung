%% Linear anatomical measurements
subplot(2, 2, 1)
% volume per node at specific level
plot(branch_volume, 'ko-', 'linewidth', 2)
xlim([1 15])
set(gca, 'XTick', 1:15)
ylabel('volume (ml)')
xlabel('levels')
title('Volume per node at the specific level')
set(gca, 'fontsize', 14)

subplot(2, 2, 2)
% volume per level
plot(branch_volume.*p.nodes_pergen, 'ko-', 'linewidth', 2)
text(2, 0.3, ['Network volume = ' num2str(sum(branch_volume.*p.nodes_pergen)) ' ml'], 'fontsize', 13)
xlim([1 15])
set(gca, 'XTick', 1:15)
ylabel('volume (ml)')
xlabel('levels')
title(' Volume per level')
set(gca, 'fontsize', 14)

subplot(2, 2, 3)
% radius
plot(branch_diameter./2, 'ko-', 'linewidth', 2)
xlim([1 15])
set(gca, 'XTick', 1:15)
ylabel('radius (cm)')
xlabel('levels')
set(gca, 'fontsize', 14)

subplot(2, 2, 4)
% length
plot(branch_length, 'ko-', 'linewidth', 2)
xlim([1 15])
set(gca, 'XTick', 1:15)
set(gca, 'YTick', 0:0.1:1)
ylabel('length (cm)')
xlabel('levels')
set(gca, 'fontsize', 14)
set(gcf, 'position', [ 1001         569        1049         770])

%% semilogy plot with anatomical measurements
levels=1:15;
subplot(2, 2, 1)
% volume per node at specific level

semilogy(levels, branch_volume, 'ko-', 'linewidth', 2)
hold on
[params,S] = polyfit(levels, log(branch_volume), 1);
slope = params(1);
intercept = params(2);
exp_func = @(data) exp(slope.*data + intercept);
semilogy(levels, exp_func(levels), '--', 'color', 'r', 'linewidth', 2)
hold off
text( 4, 0.003, ['m = ' num2str(abs(params(1)))], 'fontsize', 14)
xlim([1 15])
set(gca, 'XTick', 1:15)
ylabel('volume (ml)')
xlabel('levels')
title('Volume per node at the specific level')
set(gca, 'fontsize', 14)

subplot(2, 2, 2)
% volume per level
semilogy(levels, branch_volume.*p.nodes_pergen, 'ko-', 'linewidth', 2)
hold on
[params,S] = polyfit(levels, log(branch_volume.*p.nodes_pergen), 1);
slope = params(1);
intercept = params(2);
exp_func = @(data) exp(slope.*data + intercept);
semilogy(levels, exp_func(levels), '--', 'color', 'r', 'linewidth', 2)
hold off
text( 5, 0.02, ['m = ' num2str(abs(params(1)))], 'fontsize', 14)
text(2, 0.3, ['Network volume = ' num2str(sum(branch_volume.*p.nodes_pergen)) ' ml'], 'fontsize', 13)
xlim([1 15])
set(gca, 'XTick', 1:15)
ylabel('volume (ml)')
xlabel('levels')
title(' Volume per level')
set(gca, 'fontsize', 14)

subplot(2, 2, 3)
% radius
semilogy(levels, (branch_diameter./2), 'ko-', 'linewidth', 2)
hold on
[params,S] = polyfit(levels, log((branch_diameter./2)), 1);
slope = params(1);
intercept = params(2);
exp_func = @(data) exp(slope.*data + intercept);
semilogy(levels, exp_func(levels), '--', 'color', 'r', 'linewidth', 2)
hold off
text( 8, 0.02, ['m = ' num2str(abs(params(1)))], 'fontsize', 14)
xlim([1 15])
set(gca, 'XTick', 1:15)
ylabel('radius (cm)')
xlabel('levels')
set(gca, 'fontsize', 14)

subplot(2, 2, 4)
% length
semilogy(levels, branch_length, 'ko-', 'linewidth', 2)
hold on
[params,S] = polyfit(levels, log(branch_length), 1);
slope = params(1);
intercept = params(2);
exp_func = @(data) exp(slope.*data + intercept);
semilogy(levels, exp_func(levels), '--', 'color', 'r', 'linewidth', 2)
hold off
text( 8, 0.15, ['m = ' num2str(abs(params(1)))], 'fontsize', 14)
xlim([1 15])
set(gca, 'XTick', 1:15)
set(gca, 'YTick', 0:0.1:1)
ylabel('length (cm)')
xlabel('levels')
set(gca, 'fontsize', 14)
set(gcf, 'position', [ 1001         569        1049         770])

% %% log log
% levels=1:15;
% subplot(2, 2, 1)
% % volume per node at specific level
% 
% loglog(levels, branch_volume, 'ko-', 'linewidth', 2)
% hold on
% [params,S] = polyfit(levels, log(branch_volume), 1);
% slope = params(1);
% intercept = params(2);
% exp_func = @(data) exp(slope.*data + intercept);
% loglog(levels, exp_func(levels), '--', 'color', 'r', 'linewidth', 2)
% hold off
% text( 4, 0.003, ['m = ' num2str(abs(params(1)))], 'fontsize', 14)
% xlim([1 15])
% set(gca, 'XTick', 1:15)
% ylabel('volume (ml)')
% xlabel('levels')
% title('Volume per node at the specific level')
% set(gca, 'fontsize', 14)
% 
% subplot(2, 2, 2)
% % volume per level
% loglog(levels, branch_volume.*p.nodes_pergen, 'ko-', 'linewidth', 2)
% hold on
% [params,S] = polyfit(levels, log(branch_volume.*p.nodes_pergen), 1);
% slope = params(1);
% intercept = params(2);
% exp_func = @(data) exp(slope.*data + intercept);
% loglog(levels, exp_func(levels), '--', 'color', 'r', 'linewidth', 2)
% hold off
% text( 5, 0.02, ['m = ' num2str(abs(params(1)))], 'fontsize', 14)
% text(2, 0.3, ['Network volume = ' num2str(sum(branch_volume.*p.nodes_pergen)) ' ml'], 'fontsize', 13)
% xlim([1 15])
% set(gca, 'XTick', 1:15)
% ylabel('volume (ml)')
% xlabel('levels')
% title(' Volume per level')
% set(gca, 'fontsize', 14)
% 
% subplot(2, 2, 3)
% % radius
% loglog(levels, (branch_diameter./2), 'ko-', 'linewidth', 2)
% hold on
% [params,S] = polyfit(levels, log((branch_diameter./2)), 1);
% slope = params(1);
% intercept = params(2);
% exp_func = @(data) exp(slope.*data + intercept);
% loglog(levels, exp_func(levels), '--', 'color', 'r', 'linewidth', 2)
% hold off
% text( 8, 0.02, ['m = ' num2str(abs(params(1)))], 'fontsize', 14)
% xlim([1 15])
% set(gca, 'XTick', 1:15)
% ylabel('radius (cm)')
% xlabel('levels')
% set(gca, 'fontsize', 14)
% 
% subplot(2, 2, 4)
% % length
% loglog(levels, branch_length, 'ko-', 'linewidth', 2)
% hold on
% [params,S] = polyfit(levels, log(branch_length), 1);
% slope = params(1);
% intercept = params(2);
% exp_func = @(data) exp(slope.*data + intercept);
% loglog(levels, exp_func(levels), '--', 'color', 'r', 'linewidth', 2)
% hold off
% text( 8, 0.15, ['m = ' num2str(abs(params(1)))], 'fontsize', 14)
% xlim([1 15])
% set(gca, 'XTick', 1:15)
% set(gca, 'YTick', 0:0.1:1)
% ylabel('length (cm)')
% xlabel('levels')
% set(gca, 'fontsize', 14)
% set(gcf, 'position', [ 1001         569        1049         770])
%%

BS = res(end, 1:p.NP);
BR = res(end, p.NP+1:2*p.NP);
Btot = BS + BR;
Btot = Btot.*branch_volume';

semilogy(Btot, 'ko-', 'linewidth', 2)
xlim([1 15])
set(gca, 'XTick', 1:15)
%ylabel('Final bacterial density (CFU/ml)')
ylabel('Bacteria numbers (CFU)')
xlabel('levels')
set(gca, 'fontsize', 14)