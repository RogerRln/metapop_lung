%% Function to plot bacterial dynamics of the Immune+/Phage+ case
% we focus on showing the times when infection clears from network nodes

function plot_meatpopDyn_zoomIposPpos(res, time, p, title_label, cbar_label, y_label, x_label)
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

cm = parula; % blue-to-yellow gradient
cm = [1 1 1; cm(2:end,:)]; % white for lowest value
cmap = colormap(cm);


%indx_tinset = find(time_vec >= 31 & time_vec <= 38);
indx_tinset = find(time_vec >= 26 & time_vec <= 35);
time_inset = time_vec(indx_tinset);
imagesc(color_ind(indx_tinset,:)', [0 10]);
colormap(cmap);
cbar = colorbar;
caxis([0 10]);
ylabel(cbar, cbar_label, 'interpreter', 'latex', 'fontsize', 17)
xlabel(x_label, 'interpreter', 'latex')
ylabel(y_label, 'interpreter', 'latex')
title(title_label, 'interpreter', 'latex')
yticks(1:2:15)
ylim([0.5 15.5])
%xticks(1:12:length(time_vec(indx_tinset)))
xticks(1:20:length(time_inset))
new_tvec = round(time_inset);
%xtick = num2str(new_tvec(1:12:length(time_vec(indx_tinset))));
xtick = num2str(new_tvec(1:20:length(time_inset)));
xticklabels(xtick)
set(gca, 'fontsize', 17, 'linewidth', 2)



end