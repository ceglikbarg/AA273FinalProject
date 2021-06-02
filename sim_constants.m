%% Constants for Simulink model


% Simulation constants
constants.dt = 0.1;
constants.tend= 3600;

% System constants
constants.mu_max = 0.806/3600; % max growth rate, 1/sec
constants.Ki = 87.4; % inhibition coefficient, g/L
constants.Ks = 0.68; % saturation coefficient, g/L
constants.Yxs = 6*11.8/180; % yield coeff, g/g (11.8 g per mole Carbon * (6/180 moles Carbon per g glucose)
constants.m = (0.103/1000/3600)*180; % maintenance coeff, (mol/g/s)*g/mol = g/g/s
constants.V = 2; % volume, L
constants.C_heat = 4.314e3*constants.V; % specific heat, J/K
constants.N_power = 2*10^-10;
constants.Kq = 470000; % bacteria heat coefficient, J/(mol O2)
constants.Pp = 20; % peltier coefficient, W/A
constants.Rp = 2; % peltier resistance, ohms
constants.Np = 10; % number of peltier units at this spec
constants.hA = 1000*(100 * 0.01^2); % 1000 W/m^2/K * 10cm^2 of surface area
constants.mf_in = 0.21;
constants.S0 = 25; % concentration of substrate in feed, g/L
constants.pr = 1.2e5; % pressure in reactor, Pa

% Initial Conditions
S = 4; % g/l, initial substrate concentration
C = 2; % g/l, initial cell concentration
QO2 = 0.0198/3600; % specific O2 consumption rate for e. coli and glusoce (in mol/O2 per sec per gram)
Tc = 273.15+35; % culture temp initially, degK
Tinf = 290; % environment temp
constants.x0 = [S; C; QO2; Tc; Tinf];

% EKF Constants
constants.Q = 1E-3*eye(length(constants.x0));
constants.Q(1:2,1:2) = 1E-6*constants.dt*eye(2);
constants.Q(3,3) = 1E-15*constants.dt;
constants.R = eye(4);
constants.R(1,1) = (0.001/2)^2; % mass fraction sensor has 0.1% percent by vol error
constants.R(2,2) = (0.002/2)^2; % D sensor has 0.2% by vol error
constants.R(3,3) = (0.2/2)^2; % internal thermistor is +/- 0.2deg
constants.R(4,4) = (1/2)^2; % external thermistor is +/- 1deg

constants.mu0 = constants.x0 + [-0.0004; -0.0001; -0.0000; 0.0059; -0.0048];
constants.cov0 = 0.1*eye(length(constants.mu0));

% Targetting
constants.w = 450; % agitator speed, RPM
constants.Fo = 0.1/22.4; % moles per sec of air input (SLPS/22.4 = moles/sec)
constants.Tc_target = convtemp(98.6,'F','K'); % target temperature in K
constants.S_target = 5; % target substrate concentration, g/L
constants.u0 = [constants.w; 0; 0; constants.Fo]; % initial control input
constants.Ip_max = 10; % max amperage through peltier
constants.Fs_max = 0.5/22.4; % max flow through oxygen input