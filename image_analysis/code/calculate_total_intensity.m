%% Function to calculate total intensity signal for the top and bottom
% compartments of mice images of different immmune and phage conditions
function [intensity_top, intensity_bottom] = calculate_total_intensity(InegPneg_data, time_labels)


    intensity_top = [];
    intensity_bottom = [];

    for d = 1:length(InegPneg_data)

        top_data = InegPneg_data(d).top;
        bottom_data = InegPneg_data(d).bottom;
        pixel_intensity_top = NaN(length(top_data)-1, length(time_labels));
        pixel_intensity_bottom = NaN(length(bottom_data)-1, length(time_labels));

        for n = 2:length(top_data)

            field_names = fieldnames(top_data);

            for f = 1:length(field_names)

                indx_time = find(strcmp(field_names{f}, time_labels));
                pixel_vals_top = top_data(n).(field_names{f});
                pixel_vals_bottom = bottom_data(n).(field_names{f});

                if ~isempty(pixel_vals_top)

                    pixel_intensity_top(n-1, indx_time) = sum(sum(pixel_vals_top));
                    pixel_intensity_bottom(n-1, indx_time) = sum(sum(pixel_vals_bottom));

                end

            end

        end

        intensity_top = [intensity_top; pixel_intensity_top];
        intensity_bottom = [intensity_bottom; pixel_intensity_bottom];

    end
    
end