%% Code to plot the population dynamics of bacteria, phage and neutrophils at the node level
function figures_metapop(t,y,p)

  % Plot densities
    figure(1);
    nrows = round(p.NP/3);
    count = 0;
    for i = 1:p.NP
        subplot(nrows, 3, i)
        indx1 = [i p.NP+i 2*p.NP+i 3*p.NP+i];
        semilogy(t, y(:, indx1), 'linewidth', 3)
        ylim([1/p.branch_volume(i) 1e12])
        set(gca, 'Ytick', [1e5 1e7 1e9 1e11]);
        xlabel('Time (h)')
        ylabel('density (ml^{-1})')
        title(['Node ' int2str(i)], 'fontsize', 15)
        if ~mod(i,3)
            pos = get(gca, 'position');
            lg = legend('B_S', 'B_R', 'P', 'I', 'position', [0.9222   pos(2)    0.0590    0.0971]);
            legend boxoff
            count = count + 1;
        end
        set(gca, 'YMinorTick','off', 'linewidth', 2, 'fontweight', 'bold', 'fontsize', 13);
        xlim([0 p.T])

    end
    set(gcf, 'position', [458     1   826   827])
    % file_name = '/Users/rrodriguez77/Dropbox (GaTech)/Phage-Immune_host-pathogen project/draft_metapopulation/version2_Sep132022/figures/model_volume_theory/pop_dyn_densities.eps';
    % exportgraphics(gcf,file_name);
   
    
    % Plotting total numbers of bacteria, phage and neutrophils
    figure(2);
    nrows = round(p.NP/3);
    count = 0;
    generations = 0:p.NP-1;
    nodes_pergen = 2.^generations;
    for i = 1:p.NP
        subplot(nrows, 3, i)
        indx1 = [i p.NP+i 2*p.NP+i 3*p.NP+i];
        semilogy(t, y(:, indx1).*p.branch_volume(i), 'linewidth', 3)
        xlabel('Time (h)')
        ylabel('number of species')
        title(['Node ' int2str(i)], 'fontsize', 15)
        if ~mod(i,3)
            pos = get(gca, 'position');
            lg = legend('B_S', 'B_R', 'P', 'I', 'position', [0.9222   pos(2)    0.0590    0.0971]);
            legend boxoff
            count = count + 1;
        end
        set(gca, 'YMinorTick','off', 'linewidth', 2, 'fontweight', 'bold', 'fontsize', 13);
        ylim([1 1e11])
        yticks([1 1e3 1e6 1e9 1e11])
        xlim([0 p.T])

    end
    set(gcf, 'position', [424   206   772   908])
    % file_name = '/Users/rrodriguez77/Dropbox (GaTech)/Phage-Immune_host-pathogen project/draft_metapopulation/version2_Sep132022/figures/model_volume_theory/pop_dyn_numbers.eps';
    % exportgraphics(gcf,file_name);

end