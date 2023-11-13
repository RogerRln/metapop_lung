%% Function to plot the phage dynamics across time across the network nodes

function plot_phageDyn(res, time, p, title_label, cbar_label, y_label, x_label)

y = res;
numel = size(y,1);
mat_phagedensity = zeros(numel, p.NP);
% obtain B_tot at each time step
for i = 1:p.NP
    indx1 = 2*p.NP+i;
        for j = 1:numel
            num_phage = y(j, indx1)*p.branch_volume(i);
            if num_phage < 1
                tot_dens = 1;
            else 
                tot_dens = y(j, indx1);
            end
            mat_phagedensity(j,i) = log10(tot_dens);
        end
end

color_ind = mat_phagedensity;
new_time = time;


t_end = 50;
if new_time(end) < t_end
    time_vec =  new_time;
else
    indx_tend = find(new_time >= t_end);
    time_vec = new_time(1:indx_tend(1));
    color_ind =  color_ind(1:indx_tend(1), :);   
end

% when time include decimals, 1) round time and 2) find indexes for integers
% until t_end
round_time = round(time_vec);
[c, indx_uniq, ic] = unique(round_time, 'stable');
new_t = round_time(indx_uniq);
new_color_ind = color_ind(indx_uniq,:);

cm = parula; % blue-to-yellow gradient
cm = [1 1 1; cm(2:end,:)]; % whie for lowest value
cmap = colormap(cm);

imagesc(new_color_ind');
colormap(cmap);
cbar = colorbar;
caxis([0 12]);
ylabel(cbar, cbar_label, 'interpreter', 'latex', 'fontsize', 17)
xlabel(x_label,  'interpreter', 'latex')
ylabel(y_label,  'interpreter', 'latex')
yticks(1:2:15)
ylim([0.5 15.5])
xticks(1:5:length(new_t))
xtick = num2str(new_t(1:5:length(new_t)));
xticklabels(xtick)
title(title_label, 'interpreter', 'latex', 'fontsize', 19)
set(gca, 'fontsize', 17, 'linewidth', 2)

end