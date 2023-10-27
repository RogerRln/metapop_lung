function [chain_net, Adj_ghost] = Generate_chain_network(Ng)
%function [chain_net, Adj_ghost] = Generate_chain_network(Ng)
% Construct a Chain-metapopulation network
% Connection between patches 1-2-3-...-Ng, where Ng
% is the number of generations and each index is a generation.

    Ng = Ng;
    Adjacency_chain = zeros(Ng, Ng);

    edge_list = [];
    for i = 1:Ng
        if i < Ng
            edge_list = [edge_list; i i+1];
        end
    end

    [nr, nc] = size(edge_list);
    for i = 1:nr
        s_node = edge_list(i,1);
        t_node = edge_list(i,2);
        Adjacency_chain(s_node, t_node) =  1;
        Adjacency_chain(t_node, s_node) =  1;
        if (i > 1)
            Adjacency_chain(i, i) =  1;
        end
    end
    Adjacency_chain(end, end) =  1;
    chain_net = Adjacency_chain;
    
    
    Adj_ghost = zeros(Ng, Ng);
    for i = 1:Ng-1
       s_node = edge_list(i,1);
       t_node = edge_list(i,2);
       Adj_ghost(s_node, t_node) =  1;
    end

end



