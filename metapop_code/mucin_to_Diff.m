function [bac_diff, phage_d, n_speed] = mucin_to_Diff(mu)
%% calculate bacteria diffusion constant (DB) as a function of mucin concentration
mucin = [2.5 8];
bact_speed = [17 1]; % units um/sec
query = 2.5:0.5:8;
interp_speed = interp1(mucin, bact_speed, query);

[par_bac,S] = polyfit(query,interp_speed,1);


speed = par_bac(1)*mu + par_bac(2);
alpha = 1/3;
bac_diff = speed^2/(3*(1-alpha)); % units um^2/sec

%% calculate phage diffusion (DP) constant as a function of mucin concentration


% Low mucin concentration parameters
a_low = exp(-1.12);
b_low = 0.65;

% high mucin concentration parameters
a_high = exp(-1.13);
b_high = 1.59;

visc_W = 0.69; % mPa*s
kT = 4.28; % pN*nm
Rp = 90*1e-3; % phage radius, from nm to um
if mu <= 1
    pred_visc = a_low*mu^b_low + visc_W;
    phage_d = kT/(6*pi*pred_visc*Rp); % units um^2/sec
else
    pred_visc = a_high*mu^b_high + visc_W;
    phage_d = kT/(6*pi*pred_visc*Rp); % units um^2/sec
end


%% Calculate neutrophil speed as a function of mucin concentration

mucin = [1.5 2.5 6.5];
neutro_speed = [5 3 2.3]; % units um/min
query = 1.5:0.5:6.5;
interp_speed = interp1(mucin, neutro_speed, query);

query_high = 2.5:0.5:6.5;
interp_speed_high = interp1(mucin, neutro_speed, query_high);
[par_neutro_high,S] = polyfit(query_high,interp_speed_high,1);

query_low = 1.5:0.5:2.5;
interp_speed_low = interp1(mucin, neutro_speed, query_low);
[par_neutro_low,S] = polyfit(query_low,interp_speed_low,1);


if mu < 2.5
    n_speed = par_neutro_low(1)*mu + par_neutro_low(2);
else
     n_speed = par_neutro_high(1)*mu + par_neutro_high(2);
end

end
