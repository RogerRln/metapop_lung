% Function to calculate the neutrophil hopping rate for the influx part of
% the equation. If the function is active, neutrophils are allowed to chemotact
% at rate that is a function of neutrophil diffusion and local bacterial
% conocentration

function Tau_n = neutrophil_hopping(speed_neutro, branch_length, gamma, alpha, Btot, adj_metapop, p)
    
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
        bact_source = Btot(source(i));%*(1/branch_volume(source(i)));
        bact_target = Btot(target(i));%*(1/branch_volume(source(i)));
        Tau_n(source(i), target(i)) = gamma*((8*Diff_run_tumble)/branch_length(source(i))^2)*(1/(1 + ((bact_source/p.KB)^p.nexp)))*prob_link(source(i)) + (1-gamma)*(((8*Diff_run_tumble)/branch_length(source(i))^2)*prob_link(source(i))); % Chemo based on local bacterial concentration   
    end
  
   Tau_n = Tau_n';
   
   ind_inf = isinf(Tau_n);
   Tau_n(ind_inf) = 0;
   
end
   