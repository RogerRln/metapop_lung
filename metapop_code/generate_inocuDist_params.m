%% Code to generate the phage dose and bacterial inoculum distribution parameters
% Using a truncated gaussian function generate the parameters mean_node
% and sigma that determine the distribution of the phage dose and bacterial
% inoculum among network nodes.

%% Estas son las distribuciones (phage dose and bacterial inoculum) que uso en las simulaciones!!

load struct_variable.mat
p = struct_variable;
sigma_vec = [0 1 2 3 7 11];
initial_node = 1;
base_cdf = 0.5:1:15.5;
params_dist = [];
succ = 0;
rng(200);
for i = 1:numel(sigma_vec)
    
    node = initial_node;
    sig = sigma_vec(i);
    while node <= 15

        mean_node = node;
        mean_mu = mean_node;
        
        pd = makedist('Normal', 'mu', mean_mu , 'sigma', sig);
        t = truncate(pd, 0.5, 15.5);
        prob = zeros( 15, 1);
        
        for c = 1:numel(base_cdf)-1
            prob(c) = t.cdf(base_cdf(c+1)) - t.cdf(base_cdf(c)) ;
        end

        num_bac = prob.*1e6;
        dens = num_bac./(p.nodes_pergen.*p.branch_volume);

        if sum(dens > p.Kc/5) == 0
            params_dist = [params_dist; mean_mu sig];
            succ = succ + 1;
        end
        
        if sig == 0
            node =  node + 1;
        elseif sig > 0 && sig < 7
            node =  node + 2*sig;
        else
            node =  node + sig;
        end
        
    end
    
    if sig >= 3
        initial_node = randperm(15, 1);
        disp(initial_node)
    else
        initial_node = 1;
    end
    
end
size(params_dist)

succ = 0;
tiledlayout(6,5, 'TileSpacing', 'Compact', 'Padding', 'Compact')
prob_flag = 1;
for i = 1:size(params_dist,1)
    
    pd = makedist('Normal', 'mu', params_dist(i,1) , 'sigma', params_dist(i,2) );
    t = truncate(pd, 0.5, 15.5);
    prob = zeros( 15, 1);

    for c = 1:numel(base_cdf)-1
        prob(c) = t.cdf(base_cdf(c+1)) - t.cdf(base_cdf(c)) ;
    end
    
    num_bac = prob.*1e6;
    dens = num_bac./(p.nodes_pergen.*p.branch_volume);
    
    nexttile
    if prob_flag
        bar(prob)
        %ylim([0 1])
        ylabel('Prob.')
        xlabel('N')
        set(gca, 'fontsize', 12)
        xticks(1:2:15)
    else
        bar(dens)
        set(gca, 'YScale', 'log')
        ylabel('Density (ml$^{-1})$', 'interpreter', 'latex')
        xlabel('N')
        ylim([1e0 1e10])
        yticks([1e0 1e5  1e10])
        set(gca, 'fontsize', 12)
        xticks(1:2:15)
    end
    title(['mean = ' num2str(params_dist(i,1)) ', std = ' num2str(params_dist(i,2))])
    
    if sum(dens > p.Kc/5) == 0
       
        succ = succ + 1;
    end
    
end
set(gcf, 'position', [1001         377        1044         960])
% filename = '/Users/rrodriguez77/Dropbox (GaTech)/Phage-Immune_host-pathogen project/draft_metapopulation/version2_Sep132022/figures/variation_initial_conditions/prob_inoculum_dist.eps';
% exportgraphics(gcf, filename);


rng(24);
inocuDist_params = [];
for i = 1:size(params_dist,1)
    
    phage_params = randperm(size(params_dist,1), 3);
    inocuDist_params = [inocuDist_params; params_dist(i,:) params_dist(phage_params(1), :);...
        params_dist(i,:) params_dist(phage_params(2), :);...
        params_dist(i,:) params_dist(phage_params(3), :)];
        
    
end

% save inocuDist_params.mat inocuDist_params