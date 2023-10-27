clc;
clear;
close all

mainFolder = pwd;
folder_struct = dir(mainFolder);
subFolders = {};
count = 1;
for k = 1:length(folder_struct)
    if ((folder_struct(k).name(1) == 'N') || (folder_struct(k).name(1) == 'R')) && (folder_struct(k).name(end) == 'H')
        subFolders(count) = {folder_struct(k).name};
        count = count+1;
    end
end

mainFolder_path = mainFolder;
subfolder_struct = dir([mainFolder_path '/' subFolders{5}]);
subFolder_path = subfolder_struct(1).folder;

user_metadata = fileread([subFolder_path '/AnalyzedClickInfo.txt']);
exp_group = regexp(user_metadata, 'Comment2:\s+(.+?)\r', 'tokens');
disp(exp_group{1}{1})
mice_number = exp_group{1}{1};
mice_num = regexp(mice_number, 'B.+?\-\s(\d).(\d)', 'tokens');


% reading mice photograph to select regions of interest (ROIs)
tiff1 = Tiff([subFolder_path '/photograph.TIF'],'r');
mice_photo = read(tiff1);
mice_photo_adj = imadjust(mice_photo);
imagesc(mice_photo_adj);
title('Mice photograph')

roi_w = 89;
roi_h = 100;
roi1_pos = [98-(roi_w/2), 47, roi_w, roi_h];
roi2_pos = [182-(roi_w/2), 47, roi_w, roi_h];
roi3_pos = [269-(roi_w/2), 47, roi_w, roi_h];


roi1 = drawrectangle('Position', roi1_pos, 'color', 'r');
roi2 = drawrectangle('Position', roi2_pos, 'color', 'r');
roi3 = drawrectangle('Position', roi3_pos, 'color', 'r');
close(tiff1);

% reading luminescent image to extract pixel intensity
tiff2 = Tiff([subFolder_path '/luminescent.TIF'],'r');
mice_lumi = read(tiff2);
figure(1);
imagesc(mice_lumi);
title('Original luminescent image')
close(tiff2);

% scaling
scale = size(mice_lumi,1)/size(mice_photo,1);
border = (106-47+1)*scale;

% filtering out pixel intensity outliers
avg_intensity = mean(mean(mice_lumi));
filter = [0 1 0; 1 0 1; 0 1 0];
filtered_im = imfilter(mice_lumi, filter, avg_intensity);
filtered_im = filtered_im./sum(sum(filter));
figure(2);
imagesc(filtered_im);
title('Filtered image to find pixel outliers')


size_image = size(filtered_im,1);
index_outliers = [];
int_ratio =[];
for i = 1:size_image
    for j = 1:size_image
        
        center = [i j];
        left = [mod(i, size_image) mod(j-1, size_image)];
        right = [mod(i, size_image) mod(j+1, size_image)];
        up = [mod(i-1, size_image) mod(j, size_image)];
        below = [mod(i+1, size_image) mod(j, size_image)];
        left(left == 0) = size_image;
        right(right == 0) = size_image;
        up(up == 0) = size_image;
        below(below == 0) = size_image;
        
        neigh_intensity = [filtered_im(left(1), left(2))  filtered_im(right(1), right(2))  filtered_im(up(1), up(2))  filtered_im(below(1), below(2))];
        center_pixel = filtered_im(center(1), center(2));
        num_neigh = sum(neigh_intensity > 1.3*center_pixel);
        %itensity_ratio = neigh_intensity/center_pixel;
        if num_neigh >= 3
           index_outliers = [index_outliers;center]; 
        end
    end
end


clean_mice_lumi = mice_lumi;
for i = 1:size(index_outliers,1)
    clean_mice_lumi(index_outliers(i,1), index_outliers(i,2)) = avg_intensity;
end
figure(6);imagesc(clean_mice_lumi);

sort_pixel_intens = sort(clean_mice_lumi(:));
sort_pixel_intens = sort_pixel_intens(~(sort_pixel_intens > 1e3));
% intensity thresholds
uppp_thresh = 9.5e2;%double(mean(sort_pixel_intens(end-20:end)));
low_thresh = avg_intensity+25;

% displaying original luminescence image
figure(3);
imagesc(clean_mice_lumi, [low_thresh, uppp_thresh]);
title('Clean luminescent image')
colormap(jet);
colorbar;
roi1 = drawrectangle('Position', roi1_pos.*scale, 'color', 'r');
roi2 = drawrectangle('Position', roi2_pos.*scale, 'color', 'r');
roi3 = drawrectangle('Position', roi3_pos.*scale, 'color', 'r');



% convert to gray scale, pixel intensity range [0 1]
gray_mice_lumi = mat2gray(clean_mice_lumi, [low_thresh uppp_thresh]);
figure(4);
imagesc(gray_mice_lumi, [0, 1]);
title('Gray luminescent image')
colormap(jet);
colorbar;
roi1 = drawrectangle('Position', roi1_pos.*scale, 'color', 'r');
roi2 = drawrectangle('Position', roi2_pos.*scale, 'color', 'r');
roi3 = drawrectangle('Position', roi3_pos.*scale, 'color', 'r');


% crop gray scale image using ROIs and save

roi_vec = [roi1_pos.*scale; roi2_pos.*scale; roi3_pos.*scale];
crop_size = zeros(size(roi_vec,1), 2);
figure(5);
for n = 1:size(roi_vec, 1)
    Icrop(n).mat = imcrop(gray_mice_lumi, roi_vec(n,:));
    subplot(1, size(roi_vec,1), n)
    imagesc(Icrop(n).mat, [0 1]);
    crop_size(n,:) = size(Icrop(n).mat);
end
colormap(jet);colorbar;
sgtitle('Individual mouse (ROI sections)')

max_crop_width = max(crop_size(:, 2));
for n = 1:length(Icrop) 
   Icrop_width = size(Icrop(n).mat, 2);
   if Icrop_width < max_crop_width
       Icrop(n).mat = [Icrop(n).mat zeros(size(Icrop(n).mat,1), 1)];
   end
end

% rescaling images to be same size as the MYD88 ones and C57 black WT mice
for n = 1:length(Icrop)
   resize_img = imresize(Icrop(n).mat, 1/2);
   Icrop(n).mat = resize_img;
   subplot(1,5,n)
   imagesc(resize_img, [0 1]);
   hold on
   plot([1 size(resize_img, 2)], [border./2 border./2], '-r', 'linewidth', 1.5)
   hold off
   colormap(jet)
end

% Save data
tot_mice = str2num(mice_num{1}{2}) - str2num(mice_num{1}{1}) + 2;
mice_data_24H = zeros(size(Icrop(1).mat, 1),  size(Icrop(1).mat, 2), tot_mice);

for n = 1:tot_mice
    mice_data_24H(:,:, n) = Icrop(n).mat;
end
save mice_data_24H.mat mice_data_24H

border_24H = border./2;
save border_24H.mat border_24H