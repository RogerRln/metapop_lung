function [value, isterminal, direction] = myEventsFcn(t,y,p)
    % Locate the time when height passes through zero in a decreasing
    % direction
    % and stop integration.
    
%     value = [y(1)-1; y(16)-1];
%     isterminal = [1; 1];
%     direction = [-1; -1];

    load v_patches
    volume_patch = v_patches;   
    volume_patch = [volume_patch; volume_patch];
    
    
    NP = 15;
    CFUs = (y(1:2*NP).*volume_patch); % calculate number of bacteria
    CFUs(CFUs < 1) = 0; % find index where bacteria numbers < 1
    value = CFUs;
%     value = CFUs - 1;
%     val = CFUs - 1;
%     val(val < 0) = 0;
%     value = val;
    
    
    %value = (y(1:2*NP).*volume_patch)-1; numbers
    %value = y(1:2*NP)-1'; densities
    isterminal = ones(2*NP, 1);
    direction = ones(2*NP,1).*-1; % terminate when function is decreasing
    %direction = zeros(2*NP,1); % terminate no matter function is increasing or decreasing

end