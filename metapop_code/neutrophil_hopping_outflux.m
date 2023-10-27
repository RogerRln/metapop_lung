function Tau_n_out = neutrophil_hopping_outflux(speed_neutro, branch_length, gamma, alpha, Btot, adj_metapop, p, t)
    
    % Define your Tau matrices
    Tau_n = adj_metapop;
    
    [source, target] = find(Tau_n);
    %Btot = Btot.*branch_volume;
    
    % En el tiempo 0 cuando hay solo bacterias en el nodo 1, exp(numero de
    % bacterias target) es infinito B(target(i)) - Btot(source(i)) cuando
    % source es nodo 2 y target es nodo 1, resolver esto esperando que las
    % bacterias lleguen a los demas nodos.
    % ver chemo_neutrophils example, focal_concentration
    
    % tambien queremos usar Btot como medida de concentracion (cfu/ml) o
    % queremos usar numero total de bacteria por node, i.e.,
    % Btot*volume_node
    
    Diff_run_tumble = (speed_neutro^2)/(2*alpha);
    
    [nodes, nodes] = size(adj_metapop);
    patch_deg = repmat(4, 1, nodes); % vector with degrees of size 1 x N patches
        
    patch_deg(1) = 2;
    patch_deg(end) = 2;
    prob_link = patch_deg.^-1;
    for i = 1:length(source)
        
        bact_source = Btot(source(i)); %*(1/branch_volume(source(i)));
        bact_target = Btot(target(i)); %*(1/branch_volume(source(i)));
        %Tau_n(source(i), target(i)) = gamma*((2*speed_neutro)/branch_length(source(i)))*(1/(exp(-((bact_target-bact_source)/(p.mod*p.Kc))*gamma) + 1)) + (1-gamma)*((8*Diff_run_tumble)/branch_length(source(i))^2);
        
%         Tau_n(source(i), target(i)) = gamma*((2*speed_neutro)/branch_length(source(i)))*(1/(exp(-((bact_target-bact_source)/(p.mod*p.Kc))*gamma) + 1)) + (1-gamma)*(((8*Diff_run_tumble)/branch_length(source(i))^2)*prob_link(source(i))); % Chemotaxis based on absolute difference (modulation factor)
%         Tau_n(source(i), target(i)) = gamma*((2*speed_neutro)/branch_length(source(i)))*(1/(exp(-( (bact_target - bact_source)/(bact_source + 1) )*gamma*si) + 1)) + (1-gamma)*(((8*Diff_run_tumble)/branch_length(source(i))^2)*prob_link(source(i))); % Chemotaxis based on relative difference
        Tau_n(source(i), target(i)) = gamma*((8*Diff_run_tumble)/branch_length(source(i))^2)*(1/(1 + ((bact_source/p.KB)^p.nexp)))*prob_link(source(i)) + (1-gamma)*(((8*Diff_run_tumble)/branch_length(source(i))^2)*prob_link(source(i))); % Chemo based on local bacteria concentration        
        
%         if bact_source < 1
%             Tau_n(source(i), target(i)) = gamma*((8*Diff_run_tumble)/branch_length(source(i))^2)*prob_link(source(i)) + (1-gamma)*(((8*Diff_run_tumble)/branch_length(source(i))^2)*prob_link(source(i))); % Chemo linearized form
%         else
%             Tau_n(source(i), target(i)) = gamma*((8*Diff_run_tumble)/branch_length(source(i))^2)*(-log10(bact_source/p.Kc)/log10(p.Kc))*prob_link(source(i)) + (1-gamma)*(((8*Diff_run_tumble)/branch_length(source(i))^2)*prob_link(source(i))); % Chemo linearized form
%         end


        
    end
  
   for i = 1:size(Tau_n,1)
       indx = find(adj_metapop(i,:));
%        if (i < p.NP) && (i > 1)
%             Tau_n_out(i) = sum([Tau_n(i, indx) Tau_n(i, i+1)]);
%             %Tau_n_out(i) = mean([Tau_n(i,indx) Tau_n(i,i+1)]); % Este es el original, account for the second daughter branch
%             %Tau_n_out(i) = sum(Tau_n(i,indx))/(numel(indx)+1);
%        elseif (i == 1)
%            Tau_n_out(i) = sum(Tau_n(i, indx));
%        else
%             %Tau_n_out(i) = sum(Tau_n(i,indx))/(numel(indx)); 
%             %Tau_n_out(i) =  mean(Tau_n(i,indx)); % Este es el original
%             Tau_n_out(i) =  sum(Tau_n(i, indx)); 
%        end
        if (i == 1)
            Tau_n_out(i) = Tau_n(i, indx)*2; % Tau_n(i, indx)
        elseif (i < p.NP) && (i > 1)
            Tau_n_out(i) =  sum([Tau_n(i, i-1) Tau_n(i, i) Tau_n(i, i+1) Tau_n(i, i+1)]); % Tau_n(i, i)*4; 
        else
            Tau_n_out(i) =  sum([Tau_n(i, i-1) Tau_n(i, i)]); % Tau_n(i, i)*2;
        end
       
%        if numel(indx)
%             Tau_n_out(i) = mean(Tau_n(i,indx));
%        else
%             Tau_n_out(i) = 0;
%        end
       
   end
   
   %ind_inf = isinf(Tau_n);
   %Tau_n(ind_inf) = 0;
   
   Tau_n_out = Tau_n_out';
   %disp([t Tau_n_out(1)])
   
%    fprintf(p.fileID1,'%.10f %e\n', [Tau_n_out'; Btot']);
%    fprintf(p.fileID2,'%.10f %.10f %e %.10f %e\n', [t; Tau_n(13, 14); Btot(13); Tau_n(14, 13); Btot(14)]);


   
end
   