% function to calculate the bacterial elimination time given a specific
% bacterial extinction threshold (e.g., bacteria could become extinct when they reach 1 CFU count)

function t_clear  = time_clearance_thresh(res, time, p, threshold)


y = res;
numel = size(y,1);
mat_bactdensity = zeros(numel, p.NP);
% obtain B_tot density at each time step
for i = 1:p.NP
    indx1 = [i p.NP+i];
        for j = 1:numel
            num_bacteria_BS = y(j, indx1(1))*p.branch_volume(i);
            num_bacteria_BR = y(j, indx1(2))*p.branch_volume(i);
            if num_bacteria_BS < 1 && num_bacteria_BR < 1
                tot_dens = 0;
            elseif num_bacteria_BS >= 1 && num_bacteria_BR < 1
                tot_dens = y(j, indx1(1));
            elseif num_bacteria_BS < 1 && num_bacteria_BR >= 1
                tot_dens = y(j, indx1(2));
            else
                tot_dens = sum(y(j, indx1));
            end
            mat_bactdensity(j,i) = tot_dens;
        end
end

t_clear = zeros(p.NP,1);
for i = 1:p.NP
    
    growth_dir = diff( mat_bactdensity(:, i) );
    
    find_decline = find(growth_dir < 0);
    if find_decline(1) ~= 1 && find_decline(2) ~= 2
        decline_start = find_decline(1) + 1;
    elseif find_decline(2) ~= 2
        decline_start = find_decline(2) + 1;
    else
        decline_start = find_decline(3) + 1;
    end
    bact_dens = mat_bactdensity(decline_start:end, i);
    new_time = time(decline_start:end);
    
    ind_bacteria = find(bact_dens <= threshold(i));

    
    if ind_bacteria(1) == 1
        
        ind_previous = find(time == new_time(ind_bacteria(1))-1); % si new_time(ind_bacteria(1)) no es numero entero puede fallar esta condicion
        query_time = time(ind_previous): 0.01 :new_time(ind_bacteria(1));
        if query_time(end) ~= new_time(ind_bacteria(1))
            query_time = [query_time new_time(ind_bacteria(1))];
        end
        interp_values = interp1([time(ind_previous) new_time(ind_bacteria(1))], [mat_bactdensity(ind_previous,i) bact_dens(ind_bacteria(1))], query_time);
        times_to_clear = query_time(interp_values < threshold(i));
        t_clear(i) = times_to_clear(1);
        
    else
        
        query_time = new_time(ind_bacteria(1)-1): 0.01 :new_time(ind_bacteria(1));
        if query_time(end) ~= new_time(ind_bacteria(1))
            query_time = [query_time new_time(ind_bacteria(1))];
        end
        interp_values = interp1([new_time(ind_bacteria(1)-1) new_time(ind_bacteria(1))], [bact_dens(ind_bacteria(1)-1) bact_dens(ind_bacteria(1))], query_time);
        times_to_clear = query_time(interp_values < threshold(i));
        t_clear(i) = times_to_clear(1);
        
    end
%     t_clear(i) = new_time(ind_bacteria(1));
    
end

end