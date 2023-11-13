function dx = ode_metapopchain(t,y0,p)

dx = zeros(size(y0));

%-----------Parameters------------
% susceptible bacteria growth rate
rs = p.rs;
% resistant bacteria growth rate
rr = p.rr;
% total bacteria carrying capacity
Kc = p.Kc;
% nonlinear adsorption rate of phage:
phi = p.phi;
% power law exponent in phage infection:
g = p.g;
% immune response killing rate parameter:
ep = p.ep;
% bacterial conc. at which immune response is half as effective:
Kd = p.Kd;
% burst size of phage:
beta = p.beta;
% decay rate of phage:
w = p.w;
% maximum growth rate of immune response:
a = p.a;
% max capacity of immune response:
Ki = p.Ki;
% conc. of bacteria at which imm resp growth rate is half its maximum:
Kn = p.Kn;
% probability of emergence of phage-resistant mutation per cell division
m = p.m;
% probability of reversible mutation probability (Phage resistant -> Phage sensitive)
m2= p.m2;
% Migration rate of BS, BR
D = p.D;
% Migration rate of P
DP = p.DP;
% Adjacency matrix of patch network
A = p.A;
% number of patches
NP = p.NP;
% patch volume
PA = p.branch_volume;
% bacteria carrying capacity given choice of Area or Volume for patch size
Kc_patch = p.Kc;
% immune carrying capacity
Ki_patch = p.Ki;
% ghost daughter branches
ghost_net = p.ghost_network;

%  neutrophil transport parameters
speed_neutro = p.speed_neutro;
gamma = p.gamma;
alpha = p.alpha;
% airway length
branch_length = p.branch_length;
% adjacency matrix of metapopulation network
adj_metapop =  p.adj_metapop;

% Tau_b
Tau_b = p.Tau_b;
% Tau_p
Tau_p = p.Tau_p;


y0 = max(y0,0);
BS = y0(1:NP);
BR = y0(NP+1:2*NP);
P = y0((2*NP+1):3*NP);
I = y0((3*NP+1):4*NP);
Btot = BS + BR;

% functions of neutrophil hopping rate (in case we allow neutrophils to
% chemotact)
Tau_n_in = neutrophil_hopping(speed_neutro, branch_length, gamma, alpha, Btot, adj_metapop, p);
Tau_n_out = neutrophil_hopping_outflux(speed_neutro, branch_length, gamma, alpha, Btot, adj_metapop, p, t);

% Change in susceptible bacterial population
dBS = ((rs*BS - rs*BS.*(Btot./Kc_patch)).*(1-m)) - (((P.^g)*phi).*BS) - (ep.*I.*BS./(1 + (Btot./Kd))) + (rr*BR.*(1-(Btot./Kc_patch))*m2) - D.*BS + (Tau_b.*A)*BS;

% Change in resistant bacterial population
dBR = ((rr*BR - rr*BR.*(Btot./Kc_patch))*(1-m2)) + (rs*BS.*(1-(Btot./Kc_patch))*m) - (ep.*I.*BR./(1+(Btot./Kd))) - D.*BR + (Tau_b.*A)*BR;

% Change in phage population
dP = (((P.^g)*beta*phi).*BS) - (w*P) - DP.*P + (Tau_p.*A)*P;

% Change immune response
% immune continous models

%dI = ((a*I).*(1-(I.*(1/Ki_patch)))).*(Btot./(Btot + Kn)) - (Tau_n_out.*I) + (Tau_n_in.*p.metapop_volume)*I; % no system-level immune saturation
dI = ((a.*I).*(1-(I.*(1/Ki_patch)))).*(Btot./(Btot + Kn)).*(1-(sum(I.*PA.*p.nodes_pergen)/p.max_neutrophils)) - (Tau_n_out.*I) + (Tau_n_in.*p.metapop_volume)*I; % with system-level immune saturation
%dI = (a*I).*(Btot./(Btot + Kn)).*(1-(sum(I.*PA.*p.nodes_pergen)/p.max_neutrophils)) - (Tau_n_out.*I) + (Tau_n_in.*p.metapop_volume)*I; % no immune density saturation only system-level saturation


dx = [dBS; dBR; dP; dI];


end