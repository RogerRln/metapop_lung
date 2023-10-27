function  plot_fun_signal_intensity(intensity_data, xticks_lab, yl, title_lab, y_label, dead_alive, dead_txt_y, live_txt_y)

    if dead_alive

        cols = brewermap(2, 'Set1');
        dead = find(isnan(intensity_data(:,end)));
        live = find(~isnan(intensity_data(:,end)));

        for i = 1:length(dead)
            plot(intensity_data(dead(i),:),'linewidth', 3, 'color', cols(1,:), 'HandleVisibility', 'off')
            hold on
        end

        for i = 1:length(live)
            plot(intensity_data(live(i),:),'linewidth', 3, 'color', cols(2,:), 'HandleVisibility', 'off')
            hold on
        end

        plot(nan, nan, '-', 'linewidth',3, 'color', cols(1,:));
        plot(nan, nan, '-', 'linewidth', 3, 'color', cols(2,:));
        %plot(mean(InegPneg_pixint, 'omitnan'), 'LineWidth', 3, 'color', 'k')
        hold off
        text(1.1, dead_txt_y, ['Dead = ' num2str(length(dead))], 'fontsize', 14)
        text(1.1, live_txt_y, ['Alive = ' num2str(length(live))], 'fontsize', 14)
        legend('Dead', 'Alive', 'location', 'northwest')
        legend boxoff
        xlabel('Time (h)')
        xticks(1:length(xticks_lab))
        xticklabels(xticks_lab)
        ylabel(y_label)
        title(title_lab)
        ylim(yl)

    else
        
        if size(intensity_data,1) > 12
            
            cols = brewermap(12, 'Paired');
            extra_cols = [0 0 0; 211/255 211/255 211/255]; % black and light gray
            cols = [cols; extra_cols];
            
        else
            
            cols = brewermap(size(intensity_data,1), 'Paired');
            
        end
        
        for i = 1:size(intensity_data,1)
            plot(intensity_data(i,:),'linewidth', 3, 'color', cols(i,:))
            hold on
        end
        hold off
        %text(1.1, 30, ['Total mice = ' num2str(size(intensity_data,1))], 'fontsize', 14)
%         legend('mouse ' + string(1:size(intensity_data,1)),'location', 'northeast')
        wt_label = 'WT ' + string(1:4);
        rag2_label =  'Rag2^{-/-}Il2rg^{-/-} ' + string(5:size(intensity_data,1));
        legend([wt_label rag2_label],'location', 'eastoutside')
        legend boxoff
        xlabel('Time (h)')
        xticks(1:length(xticks_lab))
        xticklabels(xticks_lab)
        ylabel(y_label)
        title(title_lab)
        ylim(yl)

    end
    set(gca, 'fontsize', 15, 'linewidth', 1.5)
    
end