function plot_columns_clearance_IposCase(clearance_top, clearance_bottom, title_lab, ylab, p_val)


    num_top = length(clearance_top);
    num_bottom = length(clearance_bottom);

    mean_top = median(clearance_top);
    mean_bottom = median(clearance_bottom);

    std_top = std(clearance_top);
    std_bottom = std(clearance_bottom);


    min_range= 0.75;
    max_range = 1.25;
    r_top = linspace(min_range, max_range, num_top);
    r_bottom = r_top + 1;
    
    colors = brewermap(2, 'Dark2');
    
    % Rag2 mice, clearance time for Top and Bottom compartments
    for i = 5:numel(clearance_top)
       plot([1  2] , [clearance_top(i) clearance_bottom(i)], 'ko-', 'LineWidth', 1.5, 'MarkerSize', 8,'MarkerFaceColor', 'w', 'HandleVisibility', 'off')
       hold on
    end
    
     % WT mice, clearance time for Top and Bottom compartments
    for i = 1:4
       plot([1  2] , [clearance_top(i) clearance_bottom(i)], 'ko-', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'HandleVisibility', 'off')
       hold on
    end
    
    xlim([0.5 2.5])
    
    % plot mean
    plot([0.85 1.15], [mean_top mean_top], 'r-', 'LineWidth', 2.5, 'HandleVisibility', 'off')
    plot([1.85 2.15], [mean_bottom mean_bottom], 'r-', 'LineWidth', 2.5)
    
    % plot not significant line
    plot([1  2] , [75 75], '-k', 'LineWidth', 1.5, 'HandleVisibility', 'off')
    text([1.3 1.3], [77 77], ['p = ' num2str(round(p_val, 4))], 'fontsize', 15)
    
    plot([nan nan], [nan nan], 'o', 'LineWidth', 1, 'MarkerSize', 10,'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
    plot([nan nan], [nan nan], 'o', 'LineWidth', 1, 'MarkerSize', 10,'MarkerEdgeColor', 'k')
    
    hold off
    xticks([1 2])
    xticklabels({'Top', 'Bottom'})
    ylim([0 79])
    
    ylabel(ylab)
    title(title_lab)
    %text([1.3 1.3], [40 40], ['Total mice = ' num2str( max([length(clearance_top) length(clearance_bottom)] )  )], 'fontsize', 13)
    text([0.9 0.9], [-5.5 -5.5], ['N = ' num2str(num2str(num_top))], 'fontsize', 13)
    text([1.86 1.86], [-5.5 -5.5], ['N = ' num2str(num2str(num_bottom))], 'fontsize', 13)
    lg1 = legend('median', 'WT', '$Rag2^{-/-}Il2rg^{-/-}$', 'location', 'northeast');
    set(lg1,'Interpreter','latex');
    legend box off
    set(gca, 'fontsize', 15, 'box', 'off', 'LineWidth', 2)
    set(gcf, 'position', [1001         745         405         592])
    
end