function lung_volume = calculate_network_volume()
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
Ng = total_generations; % Number of generations

% Convert airway length and diameter from mm to cm, and airway volume from mm^3 to ml
branch_volume = [];
branch_length = [];
branch_diameter = [];
for i = 1:Ng
    branch_volume = [branch_volume; repmat(volume_generation(i)*1e-3, 1, 1)]; % from mm^3 to cm^3
    branch_length = [branch_length; repmat(branchlength(i)*0.1, 1, 1)]; % from mm to cm
    branch_diameter = [branch_diameter; repmat(diameter(i)*0.1, 1, 1)]; % from mm to cm
end

generations = 0:Ng-1;
nodes_pergen = 2.^generations;
nodes_pergen = nodes_pergen';

lung_volume = sum(branch_volume.*nodes_pergen); % network volume in ml

end