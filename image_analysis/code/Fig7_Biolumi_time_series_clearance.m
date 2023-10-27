clc
close all; 
clear all;

%% Load mice image data
dir_path = '../data/';
files = dir([dir_path '*.mat']);

for f = 1:length(files)
    
    load([dir_path files(f).name])
    disp(files(f).name)
end


[top_data_IposPpos_female, bottom_data_IposPpos_female] = top_bottom_fromLive(data_IposPpos_2015_02_24_female);
[top_data_IposPpos_male, bottom_data_IposPpos_male] = top_bottom_fromLive(data_IposPpos_2015_02_24_male);
[top_data_Rag2Ppos_1, bottom_data_Rag2Ppos_1] = top_bottom_fromLive(data_RagPpos_2015_02_24);
[top_data_Rag2Ppos_2, bottom_data_Rag2Ppos_2] = top_bottom_fromLive(data_RagPpos_2011_09_20_mice123);
[top_data_Rag2Ppos_3, bottom_data_Rag2Ppos_3] = top_bottom_fromLive(data_RagPpos_2011_09_20_mice456);

 
all_live_top = [top_data_IposPpos_female top_data_IposPpos_male top_data_Rag2Ppos_1 top_data_Rag2Ppos_2 top_data_Rag2Ppos_3];
all_live_bottom = [bottom_data_IposPpos_female bottom_data_IposPpos_male bottom_data_Rag2Ppos_1 bottom_data_Rag2Ppos_2 bottom_data_Rag2Ppos_3];
time_labels = fieldnames(all_live_top(1));

% Identify mice indexes that were used for calculating the infection clearance
% time
load('../data/immuno_data_top.mat')
load('../data/immuno_data_bottom.mat')
indx_alive = find(~isnan(immuno_data_bottom(:,end)));
live_top = immuno_data_top(indx_alive, :);
live_bottom = immuno_data_bottom(indx_alive, :);
times = [2 4 6 8 24 48 72];
[clearance_top, clearance_bottom, indx_significant] = find_cleartimes_test(live_top, live_bottom, 3, times);

%% Plot the bioluminescence infection signal of mice (WT + Rag2)

color_high_bottom = [255 107 0]./255;
% color_high_top = [255 215 0]./255;
color_high_top = [50 205 50]./255;
when_bar = (1:numel(indx_significant)).*7;
when_title = 1:numel(time_labels);
indx_end = numel(indx_significant)*7;
when_xticks = indx_end-6:indx_end;
figure(1)
t = tiledlayout(numel(indx_significant), numel(time_labels),'TileSpacing','Compact','Padding','Compact');
count = 1;
for i = 1:numel(indx_significant)
    
    for j = 1:numel(time_labels)
        top = all_live_top(indx_significant(i)).(time_labels{j});
        bottom = all_live_bottom(indx_significant(i)).(time_labels{j});
        nexttile
        imagesc([top;bottom], [0 1]);
        hold on
        if ismember(count, when_bar)
            cb = colorbar('fontweight', 'bold', 'fontsize', 9);
            cb.Label.String = 'Pixel intensity';
        end
        if ismember(count, [1 when_bar+1])
            ylabel(['mouse ' num2str(i)], 'fontweight', 'bold')
            yticks([1 10 20 30])
        else
            yticks([])
        end
        if ismember(count, when_xticks)
            xticks([1 10 20])
        else
            xticks([])
        end
        if ismember(count, when_title)
            t_lab = time_labels{j};
            title([t_lab(2:end) t_lab(1)])
        end
        if times(j) >= clearance_bottom(i)
            if clearance_bottom(i) ~= times(end)
                rectangle('Position', [1 size(top,1)+1 size(bottom,2)-1 size(bottom,1)-1], 'EdgeColor', color_high_bottom, 'LineWidth', 3.5)
            end
        end
        if times(j) >= clearance_top(i)
            if clearance_top(i) ~= times(end)
                rectangle('Position', [1 1 size(top,2)-1 size(top,1)-2], 'EdgeColor', color_high_top, 'LineWidth', 3.5)
            end
        end
        
        plot([1 size(top,2)], [size(top,1) size(top,1)], '--w', 'linewidth', 2)
        hold off
        axis tight;
        colormap(jet);
        count = count + 1;
        set(gca, 'fontsize', 12)
    end
    
end
xlabel(t,'Size (pixels)',  'fontsize', 13)
ylabel(t,'Size (pixels)',  'fontsize', 13)
set(gcf, 'position', [ 183          80         782        1257])
% filename = '../figures/Fig7_Biolumi_all_mice.eps';
% exportgraphics(gcf,filename)

%% Option 2: Separate WT and Rag2 data


% WT_indx = indx_significant(1:4);
% WT_top_clearance = clearance_top(1:4);
% WT_bottom_clearance = clearance_bottom(1:4);
% 
% Rag2_indx = indx_significant(5:end);
% Rag2_top_clearance = clearance_top(5:end);
% Rag2_bottom_clearance = clearance_bottom(5:end);
% 
% when_bar = (1:numel(WT_indx)).*7;
% when_title = 1:numel(time_labels);
% indx_end = numel(WT_indx)*7;
% when_xticks = indx_end-6:indx_end;
% 
% color_high_bottom = [255 107 0]./255;
% % color_high_top = [255 215 0]./255;
% color_high_top = [50 205 50]./255;
% 
% figure(2)
% %t = tiledlayout(numel(WT_indx), numel(time_labels),'TileSpacing','Compact', 'Padding','Compact');
% t = tiledlayout(numel(WT_indx), numel(time_labels),'TileSpacing','Compact');
% count = 1;
% for i = 1:numel(WT_indx)
%     
%     for j = 1:numel(time_labels)
%         %subplot(12, 7, count)
%         top = all_live_top(WT_indx(i)).(time_labels{j});
%         bottom = all_live_bottom(WT_indx(i)).(time_labels{j});
%         full_img = [top;bottom];
%         nexttile
%         imagesc([top;bottom], [0 1]);
%         hold on
%         if ismember(count, when_bar)
%             cb = colorbar('fontweight', 'bold', 'fontsize', 14);
%             cb.Label.String = 'Pixel intensity';
%         end
%         if ismember(count, [1 when_bar+1])
%             ylabel(['mouse ' num2str(i)], 'fontweight', 'bold')
%             yticks([1 10 20 30])
%         else
%             yticks([])
%         end
%         if ismember(count, when_xticks)
%             xticks([1 10 20])
%         else
%             xticks([])
%         end
%         if ismember(count, when_title)
%             t_lab = time_labels{j};
%             title([t_lab(2:end) t_lab(1)])
%         end
%         
%         if times(j) >= WT_bottom_clearance(i)
%             if WT_bottom_clearance(i) ~= times(end)
%                 rectangle('Position', [1 size(top,1)+1 size(bottom,2)-1 size(bottom,1)-1], 'EdgeColor',color_high_bottom, 'LineWidth', 3.5)
%             end
%         end
%         if times(j) >= WT_top_clearance(i)
%             if WT_top_clearance(i) ~= times(end)
%                 rectangle('Position', [1 1 size(top,2)-1 size(top,1)-2], 'EdgeColor',color_high_top, 'LineWidth', 3.5)
%             end
%         end
%         
%         plot([1 size(top,2)], [size(top,1) size(top,1)], '--w', 'linewidth', 2)
%         hold off
%         axis tight;
%         colormap(jet);
%         count = count + 1;
%         set(gca, 'fontsize', 15)
%     end
%     
% end
% xlabel(t,'Size (pixels)', 'fontsize', 15)
% ylabel(t,'Size (pixels)', 'fontsize', 15)
% set(gcf, 'position', [193         210        1001         595])
% % filename = '/Users/rrodriguez77/Dropbox (GaTech)/images_mice_NIHphage_New_desktop/IVIS_CHM2017_data_fig/figures/Biolumi_WT_mice_clearance_greenTop.eps';
% % exportgraphics(gcf,filename)
% 
% figure(3)
% when_bar = (1:numel(Rag2_indx)).*7;
% when_title = 1:numel(time_labels);
% %t = tiledlayout(numel(WT_indx), numel(time_labels),'TileSpacing','Compact', 'Padding','Compact');
% t = tiledlayout(numel(Rag2_indx), numel(time_labels),'TileSpacing','Compact');
% count = 1;
% for i = 1:numel(Rag2_indx)
%     
%     for j = 1:numel(time_labels)
%         %subplot(12, 7, count)
%         top = all_live_top(Rag2_indx(i)).(time_labels{j});
%         bottom = all_live_bottom(Rag2_indx(i)).(time_labels{j});
%         nexttile
%         imagesc([top;bottom], [0 1]);
%         if ismember(count, when_bar)
%             cb = colorbar('fontweight', 'bold', 'fontsize', 12);
%             cb.Label.String = 'Intensity';
%         end
%         if ismember(count, [1 when_bar+1])
%             ylabel(['mouse ' num2str(i)], 'fontweight', 'bold')
%         end
%         if ismember(count, when_title)
%             t_lab = time_labels{j};
%             title([t_lab(2:end) t_lab(1)])
%         end
%         hold on
%         plot([1 size(top,2)], [size(top,1) size(top,1)], '--w', 'linewidth', 1.5)
%         hold off
%         axis tight;
%         colormap(jet);
%         count = count + 1;
%         set(gca, 'fontsize', 13)
%     end
%     
% end
% xlabel(t,'Size (pixels)', 'fontsize', 13)
% ylabel(t,'Size (pixels)', 'fontsize', 13)
% set(gcf, 'position', [92          81        1075        1245])
% % filename = '/Users/rrodriguez77/Dropbox (GaTech)/images_mice_NIHphage_New_desktop/IVIS_CHM2017_data_fig/figures/Biolumi_Rag2_mice_clearance.eps';
% % exportgraphics(gcf,filename)
% 
