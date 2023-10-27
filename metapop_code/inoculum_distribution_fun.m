%% Function to generate the phage and bacterial initial inoculum given distribution parameters
%  Using the matrix (params_dist), we distribute the initial bacteria
% (tot_cfu) and phage (tot_pfu) among the network nodes

function [bac_inocu, phage_inocu] = inoculum_distribution_fun(params_dist, p, branch_volume, tot_cfu, tot_pfu)
    
    mean_bac = params_dist(1);
    sig_bac = params_dist(2);
    mean_phage =  params_dist(3);
    sig_phage = params_dist(4);
    base_cdf = 0.5:1:15.5;
    
    total_bacteria =  tot_cfu; %CFU
    total_phage = tot_pfu; %PFU
    
    % bacterial inoculum
    pd = makedist('Normal', 'mu', mean_bac , 'sigma', sig_bac);
    t = truncate(pd, 0.5, 15.5);
    
    prob = zeros( 15, 1);
    for c = 1:numel(base_cdf)-1
        prob(c) = t.cdf(base_cdf(c+1)) - t.cdf(base_cdf(c)) ;
    end

    num_bac = prob.*total_bacteria;
    bac_inocu = num_bac./(p.nodes_pergen.*branch_volume); % bacterial density CFU/ml
    
    % phage inoculum
    pd = makedist('Normal', 'mu', mean_phage , 'sigma', sig_phage);
    t = truncate(pd, 0.5, 15.5);
    
    prob = zeros( 15, 1);
    for c = 1:numel(base_cdf)-1
        prob(c) = t.cdf(base_cdf(c+1)) - t.cdf(base_cdf(c)) ;
    end

    num_phage = prob.*total_phage;
    phage_inocu = num_phage./(p.nodes_pergen.*branch_volume); % phage density PFU/ml
    
end
