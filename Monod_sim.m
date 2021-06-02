%% Simulation
clearvars; close all; clc;

%% Load constants
sim_constants;

%% Run sim
sim('Monod_Controller');

%% Plot output
plotsim(constants, tout, x_truth, mu, cov, y, u, targ_err)