%% Simulation
clearvars; close all; clc;

%% Load constants
sim_constants_delay;

%% Run sim
sim('Monod_delay_Controller');

%% Plot output
plotsim_delay(constants, tout, x_truth, mu, cov, y, u, targ_err)