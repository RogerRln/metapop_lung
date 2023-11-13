%% Function to calculate the infection clearance time using the total intensity signal of top and bottom compartments
% based on a significan intensity threshold calculate when the intensity
% signal decreases below the threshold

function [clearance_top, clearance_bottom, indx_significant] = find_cleartimes_test(top_data, bottom_data, significant_intensity, time)

    live_top = top_data;
    live_bottom = bottom_data;

    indx_significant_top = [];
    indx_significant_bottom = [];

    sig_intensity = significant_intensity;
    
    % identify mice whose infection signal in both compartments is above
    % the intensity threshold by the end of the experiment
    indx_persist_top = find(live_top(:,end) > sig_intensity);
    indx_persist_bottom = find(live_bottom(:,end) > sig_intensity);
    indx_perist_both = intersect(indx_persist_top, indx_persist_bottom);

    % Use mice with significant infection signals in any of the compartments.
    for i = 1:size(live_top,1)

        is_signif_top = sum(live_top(i, :) >= sig_intensity);
        is_signif_bottom = sum(live_bottom(i, :) >= sig_intensity);

        if is_signif_top
            indx_significant_top = [indx_significant_top; i];
        end

        if is_signif_bottom
            indx_significant_bottom = [indx_significant_bottom; i];
        end

    end
    indx_significant = [indx_significant_top; indx_significant_bottom];
    indx_significant = unique(indx_significant);

    % Use mice which have significant signals in both compartments
%     indx_significant = intersect(indx_significant_top, indx_significant_bottom);


    % filter out mice which signal never clears by the end of the
    % experiment
    indx_significant = indx_significant(~ismember(indx_significant, indx_perist_both));
    final_top = live_top(indx_significant, :);
    final_bottom = live_bottom(indx_significant, :);


    times = time;
    clear_time = [];
    for s = 1:size(final_bottom,1)

        signal = final_bottom(s, :);
        is_sig = sum(signal >= sig_intensity);
        end_intensity = signal(end);

        if is_sig == 0

            clear_time = [clear_time; 0];

        elseif end_intensity >= sig_intensity

            clear_time = [clear_time; times(end)];


        else

            peaks = find(signal >= sig_intensity);
            valleys = find(signal < sig_intensity);
            final_peak = peaks(end);
            clearance = valleys(valleys > final_peak);
            clear_t = clearance(1);
            
            query_time = times(clear_t-1): 0.1 :times(clear_t);
            interp_values = interp1([times(clear_t-1) times(clear_t)], [signal(clear_t-1) signal(clear_t)], query_time);
            times_to_clear = query_time(interp_values < sig_intensity);
            
            clear_time = [clear_time; times_to_clear(1)];


        end



    end
    clearance_bottom = clear_time;


    clear_time = [];
    for s = 1:size(final_top,1)

        signal = final_top(s, :);
        is_sig = sum(signal >= sig_intensity);
        end_intensity = signal(end);

        if is_sig == 0

            clear_time = [clear_time; 0];

        elseif end_intensity >= sig_intensity

            clear_time = [clear_time; times(end)];


        else

            peaks = find(signal >= sig_intensity);
            valleys = find(signal < sig_intensity);
            final_peak = peaks(end);
            clearance = valleys(valleys > final_peak);
            clear_t = clearance(1);
            
            query_time = times(clear_t-1): 0.1 :times(clear_t);
            interp_values = interp1([times(clear_t-1) times(clear_t)], [signal(clear_t-1) signal(clear_t)], query_time);
            times_to_clear = query_time(interp_values < sig_intensity);
            
            clear_time = [clear_time; times_to_clear(1)];


        end



    end
    clearance_top = clear_time;


end
