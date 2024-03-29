% This code generates all the model-related non-schematic figures in the
% manuscript: "Metapopulation model of phage therapy of an acute
% Pseudomonas aeruginosa lung infection", Rodriguez-Gonzalez R. et al,. (2023)

% All figures are saved in .png and .eps formats inside 'figures' folder 

% Run the script within the current PATH (not inside metapop_code nor inside
% image_analysis subfolders)

clear
clc
close all
addpath('./..')
%% Main figures

% Fig. 2 Spatiotemporal dynamics of phage therapy of a P.a. lung infection

run([pwd '/metapop_code/Fig2_metapop_chain_totalVolume.m'])
filename = [pwd '/figures/Fig2.png'];
exportgraphics(figure(1), filename, 'Resolution', 300)
filename = [pwd '/figures/Fig2.eps'];
exportgraphics(figure(1), filename)
clear
close all


% Fig. 3 Bacterial dynamics under different phage and innate immune treatments

% Fig. 3a-d
run([pwd '/metapop_code/Fig3AtoE_heatmaps_metapopDyn.m'])
filename = [pwd '/figures/Fig3a-d.png'];
exportgraphics(figure(1), filename, 'Resolution', 300)
filename = [pwd '/figures/Fig3a-d.eps'];
exportgraphics(figure(1), filename)

% Fig. 3e
filename = [pwd '/figures/Fig3e.png'];
exportgraphics(figure(3), filename, 'Resolution', 300)
filename = [pwd '/figures/Fig3e.eps'];
exportgraphics(figure(3), filename)

% Figure of phage density across the metapop network
filename = [pwd '/figures/extra_figs/Fig3phage_dens.png'];
exportgraphics(figure(2), filename, 'Resolution', 300)
filename = [pwd '/figures/extra_figs/Fig3phage_dens.eps'];
exportgraphics(figure(2), filename)

clear
close all

% Fig. 3f
run([pwd '/metapop_code/Fig3f_time_clearance_plot.m'])
filename = [pwd '/figures/Fig3f.png'];
exportgraphics(figure(1), filename, 'Resolution', 300)
filename = [pwd '/figures/Fig3f.eps'];
exportgraphics(figure(1), filename)

clear
close all

% Fig. 4 Infection dynamics as a result of varying the distribution of phage dose and bacterial inoculum

% Fig. 4a
run([pwd '/metapop_code/Fig4_plot_B0dist_vs_P0dist.m'])
filename = [pwd '/figures/Fig4a.png'];
exportgraphics(figure(1), filename, 'Resolution', 300)
filename = [pwd '/figures/Fig4a.eps'];
exportgraphics(figure(1), filename)

% Fig. 4b
filename = [pwd '/figures/Fig4b.png'];
exportgraphics(figure(3), filename, 'Resolution', 300)
filename = [pwd '/figures/Fig4b.eps'];
exportgraphics(figure(3), filename)

% Figure showing phage density across the metapop network for different
% bacterial and phage inocula distributions
filename = [pwd '/figures/extra_figs/Fig4phage_dens.png'];
exportgraphics(figure(2), filename, 'Resolution', 300)
filename = [pwd '/figures/extra_figs/Fig4phage_dens.eps'];
exportgraphics(figure(2), filename)

clear
close all

% Fig. 5 Outcomes of the robustness analysis: the probability of therapeutic success given intermediate innate immune response levels and varying phage adsorption rates

run([pwd '/metapop_code/Fig5_heatmap_adsorption_rate_nlung.m'])
filename = [pwd '/figures/Fig5.png'];
exportgraphics(figure(1), filename, 'Resolution', 300)
filename = [pwd '/figures/Fig5.eps'];
exportgraphics(figure(1), filename)
clear
close all


% Fig. 6 In vivo P.a. murine pneumonia data

run([pwd '/image_analysis/code/Fig6_Biolumi_time_series_clearance.m'])
filename = [pwd '/figures/Fig6.png'];
exportgraphics(figure(1), filename, 'Resolution', 300)
filename = [pwd '/figures/Fig6.eps'];
exportgraphics(figure(1), filename)
clear
close all

% Fig. 7 Time series of total intensity signal and the infection clearance analysis using in vivo P.a. murine pneumonia data

% Fig. 7a
run([pwd '/image_analysis/code/Fig7a_plot_tot_intensity_WT_Rag2_phage_mice.m'])
filename = [pwd '/figures/Fig7a.png'];
exportgraphics(figure(1), filename, 'Resolution', 300)
filename = [pwd '/figures/Fig7a.eps'];
exportgraphics(figure(1), filename)
clear
close all

% Fig. 7b
run([pwd '/image_analysis/code/Fg7b_plot_time_to_clearance.m'])
filename = [pwd '/figures/Fig7b.png'];
exportgraphics(figure(1), filename, 'Resolution', 300)
filename = [pwd '/figures/Fig7b.eps'];
exportgraphics(figure(1), filename)
clear
close all

%% Supplementary figures

% Fig. S1 Relationship between mucin concentration and bacterial speed and phage diffusion

% Fig. S1a
run([pwd '/metapop_code/FigS1_mucin_vs_diff.m'])
filename = [pwd '/figures/FigS1a.png'];
exportgraphics(figure(2), filename, 'Resolution', 300)
filename = [pwd '/figures/FigS1a.eps'];
exportgraphics(figure(2), filename)

% Fig. S1b
filename = [pwd '/figures/FigS1b.png'];
exportgraphics(figure(1), filename, 'Resolution', 300)
filename = [pwd '/figures/FigS1b.eps'];
exportgraphics(figure(1), filename)
clear
close all

% Fig. S2 Population dynamics at the node level under different phage and immune treatments

% Fig. S2a
run([pwd '/metapop_code/FigS2a_neutropenic_untreated.m'])
filename = [pwd '/figures/FigS2a.png'];
exportgraphics(figure(1), filename, 'Resolution', 300)
filename = [pwd '/figures/FigS2a.eps'];
exportgraphics(figure(1), filename)
clear
close all

% Fig. S2b
run([pwd '/metapop_code/FigS2b_immuno_untreated.m'])
filename = [pwd '/figures/FigS2b.png'];
exportgraphics(figure(1), filename, 'Resolution', 300)
filename = [pwd '/figures/FigS2b.eps'];
exportgraphics(figure(1), filename)
clear
close all

% Fig. S2c
run([pwd '/metapop_code/FigS2c_neutropenic_phage.m'])
filename = [pwd '/figures/FigS2c.png'];
exportgraphics(figure(1), filename, 'Resolution', 300)
filename = [pwd '/figures/FigS2c.eps'];
exportgraphics(figure(1), filename)
clear
close all

% Fig. S3 Infection clearance time difference between node i and the last node of the network as a function of the natural log of the volume ratio, V_i/V_bottom

run([pwd '/metapop_code/FigS3_volume_clearance_time_diff.m'])
filename = [pwd '/figures/FigS3.png'];
exportgraphics(figure(1), filename, 'Resolution', 300)
filename = [pwd '/figures/FigS3.eps'];
exportgraphics(figure(1), filename)
clear
close all

% Fig. S4 Robustness analysis outcome: the probability of therapeutic success given intermediate innate immune states and varying mucin levels

run([pwd '/metapop_code/FigS4_heatmap_mucin_nlung.m'])
filename = [pwd '/figures/FigS4.png'];
exportgraphics(figure(1), filename, 'Resolution', 300)
filename = [pwd '/figures/FigS4.eps'];
exportgraphics(figure(1), filename)
clear
close all

% Fig. S5 Bacterial dynamics of a well-mixed model given intermediate innate immune response levels

run([pwd '/metapop_code/FigS5_well_mixed_immunodeficiency.m'])
filename = [pwd '/figures/FigS5.png'];
exportgraphics(figure(1), filename, 'Resolution', 300)
filename = [pwd '/figures/FigS5.eps'];
exportgraphics(figure(1), filename)
clear
close all