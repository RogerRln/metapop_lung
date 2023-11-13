function [value, isterminal, direction] = myEventsFcn(t,y,p)
    % Locate the time when height passes through zero in a decreasing
    % direction
    % and stop integration.
    

    load v_patches % volume of network nodes
    volume_patch = v_patches;   
    volume_patch = [volume_patch; volume_patch];
    
    
    NP = 15;
    CFUs = (y(1:2*NP).*volume_patch); % calculate number of bacteria
    CFUs(CFUs < 1) = 0; % find index where bacteria numbers < 1
    value = CFUs;

    isterminal = ones(2*NP, 1);
    direction = ones(2*NP,1).*-1; % terminate when function is decreasing

end