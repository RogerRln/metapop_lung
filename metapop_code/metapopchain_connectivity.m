function [metapop_network, weight_mat, metapop_volume] = metapopchain_connectivity(adj_matrix, p, within_g, between_g, daupar, type)
% The metapoo_connectivity function receives as input the adjacency matrix
% that describes the metapopulation network and gives as an output a weighted matrix,
% where the weights reflect different types of dispersion behaviors, for
% example:

% 1. Homogeneous behavior: The flux between patches is proportional to the
% degree of the source patch (a.k.a. proportional to the number of connections
% of the patch where the flux is coming from).

% 2. Heterogeneous 1: The weights assigned to the adjacency matrix are
% proportional to the level of transit between patches. In this case, we
% assume the flux is greater within generation than between generation
% patches. The weights are determined by input parameters: within_g,
% between_g.

% 3. Heterogeneous 2: The weights of the edges of the metapopulation
% network have a biological component, here we consider the role of the
% mucociliary transport and so the weights connecting daughter-to-parent
% branches are larger than the weight conneting within generation or
% parent-to-daughter patches. The weights are determined by input parameters: within_g,
% between_g, and daupar. 

    if strcmp(type, 'homogeneous')  
        
        [nodes, nodes] = size(adj_matrix);
        patch_deg = repmat(4, 1, nodes); % degree 4 for middle network nodes
        patch_deg(1) = 2; % degree 2 for first network node
        patch_deg(end) = 2; % degree 2 for last network node
        
        prob_link = patch_deg.^-1; % connection probability based on node degree
        metapop_network = adj_matrix.*repmat(prob_link, p.NP, 1).*repmat(p.branch_volume', p.NP, 1)./p.branch_volume;
        
        metapop_volume = adj_matrix.*repmat(p.branch_volume', p.NP, 1)./p.branch_volume; % not considering degree only volume
        weight_mat = adj_matrix;
        
    elseif strcmp(type, 'heterogeneous') % Heterogeneous weights between links connecting different generation levels
       
        weight_mat = adj_matrix;

        gen = struct(); % Save the identifiers of nodes per generations
        [nodes, nodes] = size(weight_mat);
        for i = 1:nodes
            gen(i).members = i;
        end
             
        weight_mat = adj_matrix;
        
        % within generation weight
        wgen_w = within_g;
        % from parent-to-daughter weight
        parendau_w = between_g;
        % from daughter-to-parent weight
        daupar_w = daupar;

        for i = 1:p.NP
            
            neighbors_i = find(adj_matrix(i,:));
            generation_id = i;

            same_gen = gen(generation_id).members;
            indx_sg = ismember(neighbors_i, same_gen);
            sg = neighbors_i(indx_sg);
            weight_mat(i,sg) = wgen_w;

             % daughter-parent connections
            indx_dp = ~ismember(neighbors_i,  same_gen);
            dp = neighbors_i(indx_dp);
            
            if generation_id == 1

                indx_parent = ismember(dp, gen(1).members);
                weight_mat(i, dp(indx_parent)) = daupar_w*2; % multiply by 2 the flux that goes from daughter to parent branch

                indx_daughter = ~ismember(dp, gen(1).members);
                weight_mat(i, dp(indx_daughter)) = parendau_w; 
                
            else
                indx_parent = ismember(dp, gen(generation_id-1).members);
                weight_mat(i, dp(indx_parent)) = daupar_w*2;

                indx_daughter = ~ismember(dp, gen(generation_id-1).members);
                weight_mat(i, dp(indx_daughter)) = parendau_w;
            end


        end
        weight_sum = sum(weight_mat);
        %metapop_network = (weight_mat./repmat(weight_sum, p.NP, 1)).*repmat(p.branch_volume',p.NP,1);
        metapop_network = (weight_mat./repmat(weight_sum, p.NP, 1));

        
    end


end
