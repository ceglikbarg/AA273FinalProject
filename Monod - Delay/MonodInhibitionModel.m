%% Plot inhibited growth model for Monod inhibition
clearvars; close all; clc;

% Constants
mu_max = 0.806/3600 ; % max growth rate in 1/sec
Ki = 87.4; % inhibition constant in g/L
Ks = 0.68; % Monod constantin g/L

S = 0:0.01:750; % substrate concentrations to consider
rg = @(s) mu_max*Ki*s./(s+Ks)./(s+Ki);

figure(); hold on;
title('Monod Inhibition Model Substrate Concentration vs. Growth Rate');
xlabel('Substrate Concentration, g/L'); ylabel('Culture Cell Growth Rate, g/sec');
plot(S,rg(S));
pbaspect([1 1 1]);