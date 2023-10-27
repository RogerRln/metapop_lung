function plot_meatpopDyn_InocuDist(res, time, p, title_label)
y = res;
numel = size(y,1);
mat_bactdensity = zeros(numel, p.NP);
% obtain B_tot at each time step
for i = 1:p.NP
    indx1 = [i p.NP+i];
        for j = 1:numel
            num_bacteria_BS = y(j, indx1(1))*p.branch_volume(i);
            num_bacteria_BR = y(j, indx1(2))*p.branch_volume(i);
            if num_bacteria_BS < 1 && num_bacteria_BR < 1
                tot_dens = 1;
            elseif num_bacteria_BS >= 1 && num_bacteria_BR < 1
                tot_dens = y(j, indx1(1));
            elseif num_bacteria_BS < 1 && num_bacteria_BR >= 1
                tot_dens = y(j, indx1(2));
            else
                tot_dens = sum(y(j, indx1));
            end
            mat_bactdensity(j,i) = log10(tot_dens);
        end
end

color_ind = mat_bactdensity;
new_time = time;


t_end = 50;
if new_time(end) < t_end
    next_t = round(new_time(end));
    if next_t < new_time(end)
        next_t = next_t+1;
    end
    add_time = [next_t:t_end]';
    time_vec =  [new_time; add_time];
    new_zeros = zeros(length(add_time),p.NP);
    color_ind = [color_ind; new_zeros];
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


imagesc(new_color_ind', [0 10]);
colormap(cmap);

yticks(1:2:15)
ylim([0.5 15.5])
xticks(1:10:length(new_t))
xtick = num2str(new_t(1:10:length(new_t)));
xticklabels(xtick)
title(title_label, 'fontsize', 20, 'Interpreter', 'latex')
set(gca, 'fontsize', 17, 'linewidth', 1.5)
end