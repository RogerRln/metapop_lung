%% Compile all the image information from the experiment BalbC-RagGC-PAK_P1-2015-02-24 in a single data structure
% In the experiment Rag2−/−Il2rg−/− mice were treated with phage

clc
clear
close all

%% Load the image information from 2 to 72 hours post infection

load mice_data_2H.mat
load mice_data_4H.mat
load mice_data_6H.mat
load mice_data_8H.mat
load mice_data_24H.mat
load mice_data_48H.mat
load mice_data_72H.mat

load survival_48h.mat
load survival_72h.mat

load border_2H.mat
load border_4H.mat
load border_6H.mat
load border_8H.mat
load border_24H.mat
load border_48H.mat
load border_72H.mat

survival_48h = [1 survival_48h];
survival_72h = [1 survival_72h];


% pad zeros so all images are of the same size (n x m)
default_size = [35 24];
num_mice = size(mice_data_2H, 3);
for n = 1:num_mice
    
    if n <= size(mice_data_2H, 3)
        
        img_size = size(mice_data_2H(:,:,n));
        img = mice_data_2H(:,:,n);
        check_size = (img_size == default_size);
        time = {'H2'};
        
        if sum(check_size) == 2
            
            all_data(n).(time{1}) = img;
            
        elseif sum(check_size) == 1
            
            ind_pad = find(~check_size);
            num_to_pad = default_size(ind_pad) - img_size(ind_pad);
            if ind_pad == 1
                
                zero_pad = zeros(num_to_pad, img_size(2));
                img_wpad = [img; zero_pad];
            else
                zero_pad = zeros(img_size(1), num_to_pad);
                img_wpad = [img zero_pad];
            end
            
           all_data(n).(time{1}) = img_wpad;
            
        else
            rows_to_pad = default_size(1) - img_size(1);
            cols_to_pad = default_size(2) - img_size(2);
            zero_pad_rows = zeros(rows_to_pad, img_size(2));
            img_wpad = [img; zero_pad_rows];
            img_size_wpad = size(img_wpad);
            zero_pad_cols = zeros(img_size_wpad(1), cols_to_pad);
            img_wpad = [img_wpad zero_pad_cols];
            all_data(n).(time{1}) = img_wpad;
        end
        
        disp(size(all_data(n).(time{1})))
        top_bottom = all_data(n).(time{1});
        border = border_2H;
        border = round(border);
        if numel(border) == 1
            
            top_data(n).(time{1}) = top_bottom(1:border, :);
            bottom_data(n).(time{1}) = top_bottom(border+1:end, :);
            
        else
            
            top_data(n).(time{1}) = top_bottom(1:border(n), :);
            bottom_data(n).(time{1}) = top_bottom(border(n)+1:end, :);
            
        end
    end
    
end

default_size = [35 24];
num_mice = size(mice_data_4H, 3);
for n = 1:num_mice
    
    if n <= size(mice_data_4H, 3)
        
        img_size = size(mice_data_4H(:,:,n));
        img = mice_data_4H(:,:,n);
        check_size = (img_size == default_size);
        time = {'H4'};
        
        if sum(check_size) == 2
            
            all_data(n).(time{1}) = img;
            
        elseif sum(check_size) == 1
            
            ind_pad = find(~check_size);
            num_to_pad = default_size(ind_pad) - img_size(ind_pad);
            if ind_pad == 1
                
                zero_pad = zeros(num_to_pad, img_size(2));
                img_wpad = [img; zero_pad];
            else
                zero_pad = zeros(img_size(1), num_to_pad);
                img_wpad = [img zero_pad];
            end
            
            all_data(n).(time{1}) = img_wpad;
            
        else
            rows_to_pad = default_size(1) - img_size(1);
            cols_to_pad = default_size(2) - img_size(2);
            zero_pad_rows = zeros(rows_to_pad, img_size(2));
            img_wpad = [img; zero_pad_rows];
            img_size_wpad = size(img_wpad);
            zero_pad_cols = zeros(img_size_wpad(1), cols_to_pad);
            img_wpad = [img_wpad zero_pad_cols];
            all_data(n).(time{1}) = img_wpad;
        end
        
        disp(size(all_data(n).(time{1})))
        top_bottom = all_data(n).(time{1});
        border = border_4H;
        border = round(border);
        if numel(border) == 1
            
            top_data(n).(time{1}) = top_bottom(1:border, :);
            bottom_data(n).(time{1}) = top_bottom(border+1:end, :);
            
        else
            
            top_data(n).(time{1}) = top_bottom(1:border(n), :);
            bottom_data(n).(time{1}) = top_bottom(border(n)+1:end, :);
            
        end
    end
    
end


default_size = [35 24];
num_mice = size(mice_data_6H, 3);
for n = 1:num_mice
    
    if n <= size(mice_data_6H, 3)
        
        img_size = size(mice_data_6H(:,:,n));
        img = mice_data_6H(:,:,n);
        check_size = (img_size == default_size);
        time = {'H6'};
        
        if sum(check_size) == 2
            
            all_data(n).(time{1}) = img;
            
        elseif sum(check_size) == 1
            
            ind_pad = find(~check_size);
            num_to_pad = default_size(ind_pad) - img_size(ind_pad);
            if ind_pad == 1
                
                zero_pad = zeros(num_to_pad, img_size(2));
                img_wpad = [img; zero_pad];
            else
                zero_pad = zeros(img_size(1), num_to_pad);
                img_wpad = [img zero_pad];
            end
            
            all_data(n).(time{1}) = img_wpad;
            
        else
            rows_to_pad = default_size(1) - img_size(1);
            cols_to_pad = default_size(2) - img_size(2);
            zero_pad_rows = zeros(rows_to_pad, img_size(2));
            img_wpad = [img; zero_pad_rows];
            img_size_wpad = size(img_wpad);
            zero_pad_cols = zeros(img_size_wpad(1), cols_to_pad);
            img_wpad = [img_wpad zero_pad_cols];
            all_data(n).(time{1}) = img_wpad;
        end
        
        disp(size(all_data(n).(time{1})))
        top_bottom = all_data(n).(time{1});
        border = border_6H;
        border = round(border);
        if numel(border) == 1
            
            top_data(n).(time{1}) = top_bottom(1:border, :);
            bottom_data(n).(time{1}) = top_bottom(border+1:end, :);
            
        else
            
            top_data(n).(time{1}) = top_bottom(1:border(n), :);
            bottom_data(n).(time{1}) = top_bottom(border(n)+1:end, :);
            
        end
    end
    
end


default_size = [35 24];
num_mice = size(mice_data_8H, 3);
for n = 1:num_mice
    
    if n <= size(mice_data_8H, 3)
        
        img_size = size(mice_data_8H(:,:,n));
        img = mice_data_8H(:,:,n);
        check_size = (img_size == default_size);
        time = {'H8'};
        
        if sum(check_size) == 2
            
            all_data(n).(time{1}) = img;
            
        elseif sum(check_size) == 1
            
            ind_pad = find(~check_size);
            num_to_pad = default_size(ind_pad) - img_size(ind_pad);
            if ind_pad == 1
                
                zero_pad = zeros(num_to_pad, img_size(2));
                img_wpad = [img; zero_pad];
            else
                zero_pad = zeros(img_size(1), num_to_pad);
                img_wpad = [img zero_pad];
            end
            
            all_data(n).(time{1}) = img_wpad;
            
        else
            rows_to_pad = default_size(1) - img_size(1);
            cols_to_pad = default_size(2) - img_size(2);
            zero_pad_rows = zeros(rows_to_pad, img_size(2));
            img_wpad = [img; zero_pad_rows];
            img_size_wpad = size(img_wpad);
            zero_pad_cols = zeros(img_size_wpad(1), cols_to_pad);
            img_wpad = [img_wpad zero_pad_cols];
            all_data(n).(time{1}) = img_wpad;
        end
        
        disp(size(all_data(n).(time{1})))
        top_bottom = all_data(n).(time{1});
        border = border_8H;
        border = round(border);
        if numel(border) == 1
            
            top_data(n).(time{1}) = top_bottom(1:border, :);
            bottom_data(n).(time{1}) = top_bottom(border+1:end, :);
            
        else
            
            top_data(n).(time{1}) = top_bottom(1:border(n), :);
            bottom_data(n).(time{1}) = top_bottom(border(n)+1:end, :);
            
        end
    end
    
end

default_size = [35 24];
num_mice = size(mice_data_24H, 3);
for n = 1:num_mice
    
    if n <= size(mice_data_24H, 3)
        
        img_size = size(mice_data_24H(:,:,n));
        img = mice_data_24H(:,:,n);
        check_size = (img_size == default_size);
        time = {'H24'};
        
        if sum(check_size) == 2
            
            all_data(n).(time{1}) = img;
            
        elseif sum(check_size) == 1
            
            ind_pad = find(~check_size);
            num_to_pad = default_size(ind_pad) - img_size(ind_pad);
            if ind_pad == 1
                
                zero_pad = zeros(num_to_pad, img_size(2));
                img_wpad = [img; zero_pad];
            else
                zero_pad = zeros(img_size(1), num_to_pad);
                img_wpad = [img zero_pad];
            end
            
            all_data(n).(time{1}) = img_wpad;
            
        else
            rows_to_pad = default_size(1) - img_size(1);
            cols_to_pad = default_size(2) - img_size(2);
            zero_pad_rows = zeros(rows_to_pad, img_size(2));
            img_wpad = [img; zero_pad_rows];
            img_size_wpad = size(img_wpad);
            zero_pad_cols = zeros(img_size_wpad(1), cols_to_pad);
            img_wpad = [img_wpad zero_pad_cols];
            all_data(n).(time{1}) = img_wpad;
        end
        
        disp(size(all_data(n).(time{1})))
        top_bottom = all_data(n).(time{1});
        border = border_24H;
        border = round(border);
        if numel(border) == 1
            
            top_data(n).(time{1}) = top_bottom(1:border, :);
            bottom_data(n).(time{1}) = top_bottom(border+1:end, :);
            
        else
            
            top_data(n).(time{1}) = top_bottom(1:border(n), :);
            bottom_data(n).(time{1}) = top_bottom(border(n)+1:end, :);
            
        end
    end
    
end

exp_group_labels = 1:4;
default_size = [35 24];
num_mice = size(mice_data_48H, 3);
for n = 1:num_mice
    
    mouse_num = find(exp_group_labels == survival_48h(n)) + 1;
    if n == 1
        mouse_num = 1;
    end
    
    if n <= size(mice_data_48H, 3)
        
        img_size = size(mice_data_48H(:,:,n));
        img = mice_data_48H(:,:,n);
        check_size = (img_size == default_size);
        time = {'H48'};
        
        if sum(check_size) == 2
            
            all_data(mouse_num).(time{1}) = img;
            
        elseif sum(check_size) == 1
            
            ind_pad = find(~check_size);
            num_to_pad = default_size(ind_pad) - img_size(ind_pad);
            if ind_pad == 1
                
                zero_pad = zeros(num_to_pad, img_size(2));
                img_wpad = [img; zero_pad];
            else
                zero_pad = zeros(img_size(1), num_to_pad);
                img_wpad = [img zero_pad];
            end
            
            all_data(mouse_num).(time{1}) = img_wpad;
            
        else
            rows_to_pad = default_size(1) - img_size(1);
            cols_to_pad = default_size(2) - img_size(2);
            zero_pad_rows = zeros(rows_to_pad, img_size(2));
            img_wpad = [img; zero_pad_rows];
            img_size_wpad = size(img_wpad);
            zero_pad_cols = zeros(img_size_wpad(1), cols_to_pad);
            img_wpad = [img_wpad zero_pad_cols];
            all_data(mouse_num).(time{1}) = img_wpad;
        end
        
        disp(size(all_data(mouse_num).(time{1})))
        top_bottom = all_data(mouse_num).(time{1});
        border = border_48H;
        border = round(border);
        if numel(border) == 1
            
            top_data(mouse_num).(time{1}) = top_bottom(1:border, :);
            bottom_data(mouse_num).(time{1}) = top_bottom(border+1:end, :);
            
        else
            
            top_data(mouse_num).(time{1}) = top_bottom(1:border(n), :);
            bottom_data(mouse_num).(time{1}) = top_bottom(border(n)+1:end, :);
            
        end
    end
    
end


exp_group_labels = 1:4;
default_size = [35 24];
num_mice = size(mice_data_72H, 3);
for n = 1:num_mice
    
    mouse_num = find(exp_group_labels == survival_72h(n)) + 1;
    if n == 1
        mouse_num = 1;
    end
    
    if n <= size(mice_data_72H, 3)
        
        img_size = size(mice_data_72H(:,:,n));
        img = mice_data_72H(:,:,n);
        check_size = (img_size == default_size);
        time = {'H72'};
        
        if sum(check_size) == 2
            
            all_data(mouse_num).(time{1}) = img;
            
        elseif sum(check_size) == 1
            
            ind_pad = find(~check_size);
            num_to_pad = default_size(ind_pad) - img_size(ind_pad);
            if ind_pad == 1
                
                zero_pad = zeros(num_to_pad, img_size(2));
                img_wpad = [img; zero_pad];
            else
                zero_pad = zeros(img_size(1), num_to_pad);
                img_wpad = [img zero_pad];
            end
            
            all_data(mouse_num).(time{1}) = img_wpad;
            
        else
            rows_to_pad = default_size(1) - img_size(1);
            cols_to_pad = default_size(2) - img_size(2);
            zero_pad_rows = zeros(rows_to_pad, img_size(2));
            img_wpad = [img; zero_pad_rows];
            img_size_wpad = size(img_wpad);
            zero_pad_cols = zeros(img_size_wpad(1), cols_to_pad);
            img_wpad = [img_wpad zero_pad_cols];
            all_data(mouse_num).(time{1}) = img_wpad;
        end
        
        disp(size(all_data(mouse_num).(time{1})))
        top_bottom = all_data(mouse_num).(time{1});
        border = border_72H;
        border = round(border);
        if numel(border) == 1
            
            top_data(mouse_num).(time{1}) = top_bottom(1:border, :);
            bottom_data(mouse_num).(time{1}) = top_bottom(border+1:end, :);
            
        else
            
            top_data(mouse_num).(time{1}) = top_bottom(1:border(n), :);
            bottom_data(mouse_num).(time{1}) = top_bottom(border(n)+1:end, :);
            
        end
    end
    
end

flag_border = 1;
split(1).b = border_2H;
split(2).b = border_4H;
split(3).b = border_6H;
split(4).b = border_8H;
split(5).b = border_24H;
split(6).b = border_48H;
split(7).b = border_72H;

num_mice = length(all_data);
times = length(fieldnames(all_data));
time_labels = fieldnames(all_data);
sub_count = 1;
when_time = 1:7;
when_mouse_label = [1 8 15 22 29];
figure(1);
for n = 1:num_mice
    
    for t = 1:times
        
        if ~isempty(all_data(n).(time_labels{t}))
            
            subplot(num_mice, times, sub_count)
            imagesc(all_data(n).(time_labels{t}), [0 1])
            axis tight;
            colormap(jet);
            
            if flag_border
                if numel(split(t).b) == 1
                    hold on
                    plot([1 size(all_data(n).(time_labels{t}), 2)], [split(t).b split(t).b], '-r', 'linewidth', 1.5)
                    hold off
                else
                    
                    if strcmp(time_labels{t}, 'H72') && n > 1
                        mouse_num = n-1;
                        indx_mouse = find(survival_72h == mouse_num);
                        if numel(indx_mouse) > 1
                            indx_mouse = indx_mouse(end);
                        end
                        hold on
                        plot([1 size(all_data(n).(time_labels{t}), 2)], [split(t).b(indx_mouse) split(t).b(indx_mouse)], '-r', 'linewidth', 1.5)
                        hold off
                        
                    elseif strcmp(time_labels{t}, 'H48') && n > 1
                        mouse_num = n-1;
                        indx_mouse = find(survival_72h == mouse_num);
                        if numel(indx_mouse) > 1
                            indx_mouse = indx_mouse(end);
                        end
                        hold on
                        plot([1 size(all_data(n).(time_labels{t}), 2)], [split(t).b(indx_mouse) split(t).b(indx_mouse)], '-r', 'linewidth', 1.5)
                        hold off
                    else
                        hold on
                        plot([1 size(all_data(n).(time_labels{t}), 2)], [split(t).b(n) split(t).b(n)], '-r', 'linewidth', 1.5)
                        hold off
                    end
                end
            end
            
            if ismember(sub_count, when_time)
                lab = time_labels{t};
                title([lab(2:end) lab(1)])
            end
            if ismember(sub_count, when_mouse_label)
                if n ~=1
                    ylabel(['mouse ' num2str(n-1)],'fontweight','bold')
                else
                    ylabel('control', 'fontweight','bold')
                end
            end
            
            set(gca, 'fontsize', 15)
            sub_count = sub_count + 1;
            
        else
            sub_count = sub_count + 1;
        end
        
        
    end
end

cb = colorbar('Position', [0.92 0.1 0.02 0.83], 'fontsize', 15, 'fontweight', 'bold');
cb.Label.String = 'Pixel intensity';
set(gcf, 'position', [311          81        1556        1227])
% fig_dir =  '/Users/rrodriguez77/Dropbox (GaTech)/images_mice_NIHphage_New_desktop/IVIS_CHM2017_data_fig/figures/';
% saveas(gcf,[fig_dir 'data_RagPpos_2015_02_24'],'epsc')
%% Save and organize image data

data_RagPpos_2015_02_24.data = all_data;
data_RagPpos_2015_02_24.live= [2 3 5];
data_RagPpos_2015_02_24.dead = 4;
data_RagPpos_2015_02_24.top = top_data;
data_RagPpos_2015_02_24.bottom = bottom_data;

% save('../data/data_RagPpos_2015_02_24.mat', 'data_RagPpos_2015_02_24');


