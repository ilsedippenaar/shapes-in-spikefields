restoredefaultpath;
clear RESTOREDEFAULTPATH_EXECUTED
config

addpath external/export_fig
addpath(genpath('external/chronux_2_11'))
addpath external/hht
addpath external/kl_div
addpath external/InfoTheory
addpath external/fieldtrip
ft_defaults;

addpath plotting
addpath(genpath('analysis'))
addpath(genpath('notebooks'))
addpath(genpath('util'))
addpath tests
addpath(genpath('notebooks'))
addpath trialtable

params = [];
% Selecting
params.names = {'shape', 'saccade'};
params.length_postnoise_response = 260; % from Weiner and Ghose 2014
params.length_min_reaction_time = 200; % from Weiner and Ghose 2014
params.length_presaccade_response = 50; % this is just eyeball
% Multitaper
params.T = 0.5;
params.W = 6;
params.window_size = 0.064;
params.Fs = 1000;
% Plotting
params.freq_cutoff = 100;