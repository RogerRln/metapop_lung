% Function to calculate the neutrophil hopping rate for the outflux part of
% the equation. If the function is active, neutrophils are allowed to chemotact
% at rate that is a function of neutrophil diffusion and local bacterial
% conocentration

function Tau_n_out = neutrophil_hopping_outflux(speed_neutro, branch_length, gamma, alpha, Btot, adj_metapop, p, t)
    
    % Define your Tau matrices
    Tau_n = adj_metapop;
    
    [source, target] = find(Tau_n);
    
    Diff_run_tumble = (speed_neutro^2)/(2*alpha);
    
    [nodes, nodes] = size(adj_metapop);
    patch_deg = repmat(4, 1, nodes); % vector with degrees of size 1 x N patches
        
    patch_deg(1) = 2;
    patch_deg(end) = 2;
    prob_link = patch_deg.^-1;
    
    for i = 1:length(source)
        bact_source = Btot(source(i)); %*(1/branch_volume(source(i)));
        bact_target = Btot(target(i)); %*(1/branch_volume(source(i)));
        Tau_n(source(i), target(i)) = gamma*((8*Diff_run_tumble)/branch_length(source(i))^2)*(1/(1 + ((bact_source/p.KB)^p.nexp)))*prob_link(source(i)) + (1-gamma)*(((8*Diff_run_tumble)/branch_length(source(i))^2)*prob_link(source(i))); % Chemo based on local bacterial concentration         
    end
  
   for i = 1:size(Tau_n,1)
       indx = find(adj_metapop(i,:));

        if (i == 1)
            Tau_n_out(i) = Tau_n(i, indx)*2; % Tau_n(i, indx)
        elseif (i < p.NP) && (i > 1)
            Tau_n_out(i) =  sum([Tau_n(i, i-1) Tau_n(i, i) Tau_n(i, i+1) Tau_n(i, i+1)]); % Tau_n(i, i)*4; 
        else
            Tau_n_out(i) =  sum([Tau_n(i, i-1) Tau_n(i, i)]); % Tau_n(i, i)*2;
        end
       

       
   end
    
   Tau_n_out = Tau_n_out';

   
end
   