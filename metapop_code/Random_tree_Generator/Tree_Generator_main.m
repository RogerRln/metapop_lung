% This program generate an ensnmble of random trees, the outputs are number
% of heminode, total nodes, adjacency matrices and ....
clear all; clc;
it=1000;                 % Number of trails
p0=0.0;                 % Probability of zero branching 
m0=0; m1=1; m2=2;        % The bracings for full and genreal binary trees
ng=4;                    % The maximum number of allowed generations
nd=2;                    % Mmaximum nmber of daughters 
nl= nd^ng;               % Maximum number of perhiphrals (Cayley tree)
nn=(nd^(ng+1)-1)/(nd-1); % Maximum total number of the nodes (Cayley tree)
%binomail1=makedist('Binomial',nd,p0); % make the distribution for Binomial p0, and N=nd
%truncate_binomial=truncate(binomail1,1,inf);           % truncate(bi1,lower,upper), truncated binomial distrubtion (exlcude the zero branching). 

for ii=1:it
   %% Choose your favorite branching from the following:
 B=Full_Binary_Branching(nn,m0,m2,p0);               % Full Binary trees, Branching is 0(m0) or 2(m2), (set ng=4)
 %B=General_Binary_Branching(nn,m0,m1,m2,p0);        % General Binary trees, Branching is 0 (m0),1(m1), or 2(m2), (set ng=3)
% B=Uniform_Branching(nn,nd);                        % Uniform Branching, zero branching start after g>=2, nd=4,ng=4
 %B=Binomial_Branching(nn,nd,truncate_binomial,p0);  % Binomail branching, zero brnaching is start after g>=2, nd=4,ng=4
%% Adjacency Matrix 
 [adj,nh,S,n1,node]=adjacency_matrix_generator(B,ng);     % constructing the adjacency matrix from the branching
 %Save_adjacency2file(ii, nh,node,adj);                  % save the trees as file_number.tree
%% Laplacian 
 deg=zeros(node,node);
 for j=1:node 
     deg(j,j)=sum(adj(j,:)); %Degrees of each node
 end 
 laplacian=deg-adj;          %Laplacian matrix
 eig1=eig(laplacian);
 lap1(ii)=eig1(2);           %Save second eigenvalue of the laplaxcian 
%% Saving tree properties to arryas
   r1=node-S(ng); 
   k2=0; aa=find(B(1:r1)==0); k2=length(aa);    
   
   nh_inside(ii)=k2;          %The number of heminodes that are in g<ng 
   nh_g4(ii)=S(ng);           %The number of node in g=ng, 
   T(ii,1:ng)=S(1:ng);        %The population of shells 1-ng
   T(ii,ng+1)=nh;             %The total number of heminodes
   T(ii,ng+2)=node;           %The total number of nodes
   T(ii,ng+3)=nh/node;        %The ratio R=nh/node
   T(ii,ng+4)=nh/(node*node); %The ratio nh/(node^2)
   N1(ii)=node; H1(ii)=nh; 
end


%% Visualize symmetrical tree
clf
% hist(N1,100);  
% title('Histogram of the total number of the nodes')

g = graph(adj);
num_branches = size(g.Edges(:,1),1);

edge_labels = {};
for i = 1:num_branches
    edge_labels{i} = ['b', int2str(i)];
end

figure(1);
plt = plot(g, 'Edgelabel', edge_labels, 'linewidth', 2);
plt.MarkerSize = 8;
plt.NodeFontSize = 9;
plt.EdgeFontSize = 11;
plt.EdgeFontWeight = 'bold';


%% Generate metapopulation network from a symmetrical tree (w/ dichotomous branching pattern)

% obtain branches (Edges) from symmetrical tree
edges = g.Edges(:,1);
edges = table2array(edges);

% construct the metapopulation network
[num_nodes, ncol] = size(edges);

edges = [edges (1:num_nodes)'];% add a third column to the edge list to include the node number of that branch

% [nr, nc] = size(adj);
% edges = [];
% c = 1;
% for i = 1:nr
%     for count = 1:nc-c
%         j = count + c;
%         if (adj(i,j) && adj(j,i) == 1)
%             edges = [edges; i j];
%         end
%     end
%     c = c + 1;
%     j = 0;
% end

% Identify sister branches in the edge list, i.e., branches with the same parent node
% unq_parentnodes = unique(edges(:,1));
% sister_branches = [];
% for u = 1:length(unq_parentnodes)
%     indx = find(edges(:,1) == unq_parentnodes(u));
%     if length(indx) > 1
%         sister_branches = [sister_branches; edges(indx,2)'];
%             
%         
%     end
% end
%edges = [edges;sister_branches; 1 3];

source_nodes = edges(:,1);
target_nodes = edges(:,2);
new_edges = struct();
% Identify branch interactions and save them in new_edges struct variable
for i = 1:num_nodes
    branch = edges(i,1:2);
    indx_1 = find(source_nodes == branch(1));
    indx_2 = find(target_nodes == branch(1));
    indx_3 = find(source_nodes == branch(2));
    indx_4 = find(target_nodes == branch(2));
    
    indx_all = [indx_1; indx_2; indx_3; indx_4];
    
    uniq_indx = unique(indx_all);
    uniq_indx = uniq_indx(uniq_indx ~= i);
    new_edges(i).interactions = uniq_indx;

end

% Fill the adjacency matrix of the metapopulation network given branch
% interactions
adj_metapop = zeros(num_nodes, num_nodes);
for i = 1:num_nodes
    s_node = i;
    t_node = new_edges(i).interactions;
    adj_metapop(s_node, t_node) = 1;
end


%% Visualize metapopulation network

% New Graph with braches as nodes
% pair_nodes  = [1 2; 2 3; 2 4; 3 4; 3 5; 3 6; 5 6; 4 7; 4 8; 7 8;];
%patch_net = graph(metapop_edges(:,1), metapop_edges(:, 2));
metapop_net = graph(adj_metapop);
num_nodes = size(metapop_net.Nodes,1);
node_labels = {};
for i = 1:num_nodes
    node_labels{i} = ['b', int2str(i)];
end

figure(2);
plt = plot(metapop_net, 'NodeLabel', node_labels, 'linewidth', 2, 'layout', 'layered');
plt.MarkerSize = 12;
plt.NodeFontSize = 13;
%plt.EdgeFontSize = 11;
plt.NodeFontWeight = 'bold';

adjmat_patch = adjacency(metapop_net);
full(adjmat_patch)