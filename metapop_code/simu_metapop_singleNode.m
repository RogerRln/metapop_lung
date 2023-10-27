%% Function to simulate the well-mixed model given different immune levels
% 1. We assign values to parameters
% 2. We simulate the well-mixed model
% 3. We return the population dynamics of B_S, B_R, P and I and the time vector

% input: (1) B - bacterial inoculum (CFU), (2) P - phage dose (PFU)
%        (3) I - initial number of neutrophils (4) max_neutrophil_num - maximum allowed number of lung
%         neutrophils
%        (5) simu_time - simulation time

% output: (1)time- vector with simulated times
%         (2)res - matrix of size length(time) x 60 with
%         results of the population levels of BS, BR, P, and I
%         (3)p - structure variable with parameter values

function [time, res, p] = simu_metapop_singleNode(B, P, I, max_neutrophil_num, simu_time)
%% Read the anatomical information of mice airways from Counter et al. 2013

CT_Counter = readtable('./Counter_2012_airways.csv');
CT_Counter = CT_Counter(1:end,1:5);
CT_Counter = table2array(CT_Counter);
[r,c] = size(CT_Counter);
anatomy_airways = CT_Counter;


generations = anatomy_airways(:,1);
diameter=  anatomy_airways(:,2); % airway in mm
diameter_std = anatomy_airways(:,3);
branchlength = anatomy_airways(:,4); % airway length in mm
branchlength_std = anatomy_airways(:,5);

% The total lung capacity of mice is 1 ml (i.e., lung volume is 1ml) According to Irvin et al. 2003
airways = 2.^(generations); % number of airways per generation, assuming dichotomous branching pattern
total_lungvolume = 1; % ml

%% Calculate the number of generations of the metapopulation network

% We do so by calculating the network volume (sum of the volume
% of individual airways at the generation level) and comparing it
% to the total lung capacity

num_airways = readtable('./digitalize_number_airways/num_airways_pergen.csv');
num_airways = table2array(num_airways);
num_airways = round(num_airways);

rad = diameter./2;
volume_generation = pi*rad.^2.*branchlength; % in mm^3
lung_volume = total_lungvolume*1e3;% 1 cm^3 = 1e3 mm^3, According to Irvin et al. 2003


% Find the number of generations of the metapop network
lung_symmetrical = 0;
i = 0;
while lung_symmetrical <= lung_volume
    lung_symmetrical = (2^i)*volume_generation(i+1) + lung_symmetrical;
    if lung_symmetrical >= lung_volume
        break
    end
    i = i+1;  
end

last_gen = i-1; % Last generation giving a network volume <= total lung capacity
total_generations = last_gen + 1; % total number of nodes of the metapulation network

% Generate network structure
Ng = total_generations; % Number of generations
[adj_metapop, adj_ghost] = Generate_chain_network(Ng); % Chain metapopulation network with Ng nodes
[nodes, nodes] = size(adj_metapop);
p.NP = nodes; % number of nodes


% Convert airway length and diameter from mm to cm, and airway volume from mm^3 to ml
branch_volume = [];
branch_length = [];
branch_diameter = [];
for i = 1:Ng
    branch_volume = [branch_volume; repmat(volume_generation(i)*1e-3, 1, 1)]; % from mm^3 to cm^3
    branch_length = [branch_length; repmat(branchlength(i)*0.1, 1, 1)]; % from mm to cm
    branch_diameter = [branch_diameter; repmat(diameter(i)*0.1, 1, 1)]; % from mm to cm
end

p.branch_volume = branch_volume; % individual airway volume in ml
p.branch_length = branch_length; % branch length in cm
p.adj_metapop = adj_metapop; % network structure

%% Assign values to parameters
% We save the parameters in the structure variable 'p'

% % calculate bacteria and phage diffusion constants based on mucin level
% mucin_level = 2.5; % mucin concentration
% [bac_d, phage_d, neutro_s] = mucin_to_Diff(mucin_level);

% Calculate Dispersion of bacteria given a diffusion coefficient, t = x^2/diff
diff_bact = 0; % um^2/sec
diff_bact = diff_bact*1e-8*3600; % conversion to cm^2/h
time_inbranch = ((branch_length./2).^2)./(2*diff_bact); % hours
dispersion_bact = time_inbranch.^-1; % dispersion rate of bacteria in each patch, h^-1

% Calculate Dispersion of phage given its diffusion coefficient
diff_phage = 0; % um^2/sec
diff_phage = diff_phage*1e-8*3600; % conversion to cm^2/h
time_inbranch = ((branch_length./2).^2)./(2*diff_phage); % hours
dispersion_phage = time_inbranch.^-1; % dispersion rate of phage, h^-1

% Neutrophil speed based on mucin level
speed_neutro = 0; % um/min
speed_neutro = speed_neutro*1e-4*60; % cm/h


% Parameter values

% Convert concentration units of some parameters from Roach et al., 2017
% from grams (lung mass 0.135 g) to volume (network volume 0.9 ml)  

generations = 0:14;
nodes_pergen = 2.^generations;
nodes_pergen = nodes_pergen';
p.nodes_pergen = nodes_pergen;

lung_volume = sum(branch_volume.*p.nodes_pergen); % network volume in ml
p.lung_volume = lung_volume;
lung_mass = 0.135; % lung weight in grams

% susceptible bacteria growth rate
p.rs = 0.75;
% resistant bacteria growth rate
p.rr = 0.675;
% total bacteria carrying capacity in the largest branch (trachea)
%p.Kc = 1e10;
p.Kc = (1e10*lung_mass)/lung_volume;
% power law exponent in phage infection:
%p.g = 0.62; % 0.6
p.g = 0.6; % 0.6
% nonlinear adsorption rate of phage:
p.phi = 5.4e-8*((lung_volume/lung_mass)^p.g);
%p.phi = 5e-6;
% immune response killing rate parameter:
%p.ep = 8.2e-8;
p.ep = 5.4694e-07; % el que uso
% bacterial conc. at which immune response is half as effective:
p.Kd = (4.1e7*lung_mass)/lung_volume;
%p.Kd = p.Kc./244;
% burst size of phage:
p.beta = 100;
% decay rate of phage:
p.w = 0.07;
% maximum growth rate of immune response:
p.a = 0.97;
% max capacity of immune response
p.Ki = (2.4e7*lung_mass)/lung_volume;
% conc. of bacteria at which imm resp growth rate is half its maximum:
p.Kn = 1e7;
%p.Kn = p.Kc./1e3;
% probability of emergence of phage-resistant mutation per cell division
p.m = 2.85e-8; %
% probability of reversible mutation probability (Phage resistant -> Phage sensitive)
p.m2 = 2.85e-8;
% Migration rate of bacteria
p.D = dispersion_bact;
% Migration rate of phage
p.DP = dispersion_phage;

% neutrophil speed 
p.speed_neutro = speed_neutro;
p.gamma = 1;
p.alpha = 1;
p.KB = 1e5;
p.nexp = 1;

[Tau_b, Tau_p, Tau_n] = dispersion_calculation(diff_bact, diff_phage, branch_length, adj_metapop);
p.Tau_b = Tau_b;
p.Tau_p = Tau_p;
p.Tau_n = Tau_n;

% Select the type of dispersion that is ocurring in the metapopulation
% network:

% 1. Homogeneous dispersion: The metapopulation network is a non-weighted network.
% The fluxes between patches are proportional to the
% degree of the source patches.

% 2. Heterogeneous 1: The weights assigned to the metapop network are
% proportional to the level of transit between patches. In this case, we
% assume the transit is greater within generation than between generation
% patches. The weights are determined by input parameters: within_g,
% between_g.

% 3. Heterogeneous 2: The weights of the edges of the metapopulation
% network considers a biological component, a.k.a. the role of the
% mucociliary transport and so the weights connecting daughter to parent patches (i.e.,
% going from lower to higher generations) are larger than the weights conneting within generation or
% parent-to-daughter patches. The weights are determined by input parameters: within_g,
% between_g, and daupar. 


% within generation weight
wgen_w = 5;
% between generations weight (or parent-to-daughter when 'heterogeneous_2')
betgen_w= 3;
% daughter-to-parent weight (considering mucociliary transport)
daupar_w = 20;
% type of dispersion ('homogeneous', 'heterogeneous')
type = 'homogeneous';

% At the moment, the model is configured to work only with the homogeneous case
[metapop_network, weight_mat, metapop_volume] = metapopchain_connectivity(adj_metapop, p, wgen_w, betgen_w, daupar_w, type);
p.type = type;
p.weight_mat = weight_mat;



% The weight matrix determines the weights of links connecting focal patch
% 'i' with source patch 'j' but the TRANSPOSE of that weight matrix determines how much
% flux is coming from source patch 'j' to focal patch 'i'. It only applies
% to the 'heterogeneous' case.


nz_indexes = find(adj_ghost);
if strcmp(type, 'homogeneous')
    for i = 1:Ng
        if i < Ng
        metapop_network(i, i+1) =  metapop_network(i, i+1)*2; % multiply by two because the parent branch receives two times the flux from a single daughter branch
        metapop_volume(i, i+1) = metapop_volume(i, i+1)*2;
        end
    end
    p.A = metapop_network;
    adj_ghost(nz_indexes) = metapop_network(nz_indexes);
    p.ghost_network = adj_ghost;
    p.metapop_volume = metapop_volume;
    
elseif strcmp(type, 'heterogeneous')
    metapop_network = metapop_network';
    metapop_network = metapop_network.*repmat(p.branch_volume',p.NP,1);
    p.A = metapop_network;
    adj_ghost(nz_indexes) = p.A(nz_indexes);
    p.ghost_network = adj_ghost;
end

%% Simulate the metapopulation model of a P.a. lung infection

% Initialize state variable vectors
B_S0 = zeros(p.NP,1); % phage-sensitive bacteria
B_R0 = zeros(p.NP,1); % phage-resistant bacteria
P0 = zeros(p.NP,1);   % phage
I0 = zeros(p.NP,1);   % neutrophils


p.T = simu_time; % hours, simulation end time
p.max_neutrophils = max_neutrophil_num; % maximum number of lung neutrophils allowed

% Baseline Immune density level
I0(1) = I/lung_volume; % cells/ml

% Initial bacterial density
B_S0(1) = B/lung_volume; % cfu/ml

y0 = [B_S0; B_R0; P0; I0];

options = odeset('Events', @myEventsFcn_single_node, 'RelTol', 1e-6, 'AbsTol', 1e-9, 'MaxStep', 0.1);
tic
[t,y] = ode45(@(t,y) ode_metapopchain_wellmixed(t,y,p), [0:1:2], y0, options);

BS2 = y(end, 1:p.NP)';
BR2 = y(end, p.NP+1:2*p.NP)';
P2 = y(end, (2*p.NP+1):3*p.NP)';
I2  = y(end,(3*p.NP+1):4*p.NP)';
% add phage 2 hours after bacterial infection
P2(1) = P/lung_volume; % pfu/ml

tspan = 2:p.T;

yi = [BS2; BR2; P2; I2];

[t3, y3, te, ye, ie] = ode45(@(t,y) ode_metapopchain_wellmixed(t,y,p), tspan, yi, options);

if te
    disp([te t3(end) ie])
end


count = 0;
time_loop = [];
res_loop = [];
currentTime = t3(end);
tic
if currentTime < p.T-1 % Bacterial pop died before end of simulation
    y_pre = y3;
    while currentTime < p.T-1

        BS3 = y_pre(end, 1:p.NP)';
        for i = 1:length(BS3)
            if BS3(i)*lung_volume < 1
                BS3(i) = 0;
            end
        end
        BR3 = y_pre(end, p.NP+1:2*p.NP)';
        for i = 1:length(BR3)
            if BR3(i)*lung_volume < 1
                BR3(i) = 0;
            end
        end

        P3 = y_pre(end, (2*p.NP+1):3*p.NP)';
        I3  = y_pre(end,(3*p.NP+1):4*p.NP)';
        tspan3 = (currentTime+0):p.T;
        yii = [BS3;BR3;P3;I3];
        
        % This is the new global extinction threshold that stops
        % the simulation once the bact numbers are < 1 in all nodes
        if sum(BS3) == 0 && sum(BR3) == 0 
            time_loop = [time_loop; currentTime+.01];
            res_loop = [res_loop; yii'];
            p.T = currentTime+.01;
            break
        end

        % simulating diff eq
        [t4, y4, te, ye, ie] = ode45(@(t,y) ode_metapopchain_wellmixed(t,y,p), tspan3, yii, options);
        y_pre = y4;
        time_loop = [time_loop; t4(1:end)];
        res_loop = [res_loop; y4(1:end,:)];
        currentTime = t4(end);
        count = count +1;
        if te
            disp([te t4(end) ie])
        end
        
    end
    time = [t(1:end-1); t3; time_loop]; % change here -> t3(1:end-1)
    res = [y(1:end-1,:); y3; res_loop]; % change here -> y3(1:end-1, :)
    
else
    time = [t(1:end-1); t3];
    res = [y(1:end-1, :); y3]; 
end



BS_final = res(end, 1); % density
BR_final = res(end, p.NP+1); % density
num_bs = BS_final*lung_volume; % numbers
num_br = BR_final*lung_volume; % numbers
extinct_bs = num_bs < 1;
extinct_br = num_br < 1;
BS_final(extinct_bs) = 0;
BR_final(extinct_br) = 0;

Btot = BS_final + BR_final;
Btot = Btot*lung_volume;
disp(['Total bacteria (CFU) after ' num2str(p.T) 'h: ' num2str(Btot)])

end