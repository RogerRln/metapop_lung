%% Function to retrive image data only from alive mice (i.e., mice that survive the whole experiment)
function [top_data, bottom_data] = top_bottom_fromLive(all_data)

    data = all_data;
    live_mice =  data.live;
    top_data = [];
    bottom_data = [];
    for i = 1:numel(live_mice)
        top_data = [top_data data.top(live_mice(i))];
        bottom_data = [bottom_data data.bottom(live_mice(i))];      
    end

end
