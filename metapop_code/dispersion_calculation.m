% Function to calculate phage and bacteria hopping rates as a function of
% branch length and their diffusion constants.

function [Tau_b, Tau_p, Tau_n] = dispersion_calculation(diff_bact, diff_phage, branch_length, adj_metapop)
    
    % Define your Tau matrices
    Tau_b = adj_metapop;
    Tau_p = adj_metapop;
    Tau_n = adj_metapop;
    
    [source,target] = find(Tau_b);
    
    for i = 1:length(source)

            % use when considering crossing the whole branch length instead
            % of half the branch length
%             Tau_b(source(i), target(i)) = (branch_length(source(i))^2/diff_bact);
%             Tau_p(source(i), target(i)) = (branch_length(source(i))^2/diff_phage);
%             Tau_n(source(i), target(i)) = (branch_length(source(i))/speed_neutro);

            % use when considering cross half the branch length
            Tau_b(source(i), target(i)) = ((branch_length(source(i))/2)^2)/(2*diff_bact);
            Tau_p(source(i), target(i)) = ((branch_length(source(i))/2)^2)/(2*diff_phage);
            %Tau_n(source(i), target(i)) = (branch_length(source(i)/2)/speed_neutro);
    end
  
   Tau_b = Tau_b';
   Tau_p = Tau_p';
   
   Tau_b = Tau_b.^-1;
   Tau_p = Tau_p.^-1;
   
   ind_inf = isinf(Tau_b);
   Tau_b(ind_inf) = 0;
   Tau_p(ind_inf) = 0;
   
end
   