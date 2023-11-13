function [value, isterminal, direction] = myEventsFcn_single_node(t,y,p)
    % Locate the time when height passes through zero in a decreasing
    % direction
    % and stop integration.
    

    load v_patches % volume network nodes
    volume_patch = v_patches;
    
    nodes_pergen = 2.^(0:14);
    nodes_pergen = nodes_pergen';
    volume_patch = repmat(sum(volume_patch.*nodes_pergen), 15,1);
    
    volume_patch = [volume_patch; volume_patch];
    
    
    NP = 15;
    CFUs = (y(1:2*NP).*volume_patch); % calculate number of bacteria
    CFUs(CFUs < 1) = 0; % find index where bacteria numbers < 1
    value = CFUs;

    isterminal = ones(2*NP, 1);
    direction = ones(2*NP,1).*-1; % terminate when function is decreasing

end