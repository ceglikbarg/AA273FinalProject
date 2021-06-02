%% Monod EKF Testing
clearvars; clc; close all;
rng(1);

%% Sim time
dt = 0.1;
tend = 3600; 
tspan = 0:dt:tend;
tsteps = length(tspan);
%% Constants
mu_max = 0.806/3600; % max growth rate, 1/sec
Ki = 87.4; % inhibition coefficient, g/L
Ks = 0.68; % saturation coefficient, g/L
Yxs = 6*11.8/180; % yield coeff, g/g (11.8 g per mole Carbon * (6/180 moles Carbon per g glucose)
m = (0.103/1000/3600)*180; % maintenance coeff, (mol/g/s)*g/mol = g/g/s
V = 2; % volume, L
C_heat = 4.314e3*V; % specific heat, J/K
N_power = 2*10^-10;
Kq = 470000; % bacteria heat coefficient, J/(mol O2)
Pp = 20; % peltier coefficient, W/A
Rp = 2; % peltier resistance, ohms
Np = 10; % number of peltier units at this spec
hA = 1000*(100 * 0.01^2); % 1000 W/m^2/K * 10cm^2 of surface area
mf_in = 0.21;
S0 = 25; % concentration of substrate in feed, g/L
pr = 1.2e5; % pressure in reactor, Pa

%% Inputs
w = [475*ones(1,round(tsteps/3+1,0)), 475*ones(1,round(tsteps/3,0)), 475*ones(1,round(tsteps/3,0))].'; % rpm of agitator
Fs = 0.001*[0.5*ones(1,round(tsteps/3+1,0)), 0.1*ones(1,round(tsteps/3,0)), 0*ones(1,round(tsteps/3,0))].'; % substrate feed volumetric rate, L/s (1 cc per sec)
Ip = [-5*ones(1,round(tsteps/3+1,0)), 3*ones(1,round(tsteps/3,0)), 4*ones(1,round(tsteps/3,0))].'; % Peltier current
Fg = ones(tsteps,1); %((50/1000)*(1 + tspan./(8*3600))/60).'; %sccm to slps
u = [w, Fs, Ip, Fg];

%% Simulation for truth, measurements
QO2 = 0.0198/3600; % specific O2 consumption rate for e. coli and glusoce (in mol/O2 per sec per gram)
Tc = 273.15+37; % culture temp initially, degK
Tinf = 290; % environment temp
x0 = [4; 2; QO2; Tc; Tinf];

Q = 1E-3*eye(length(x0));
Q(1:2,1:2) = 1E-6*dt*eye(2);
Q(3,3) = 1E-15*dt;
R = eye(4);
R(1,1) = (0.001/2)^2; % mass fraction sensor has 0.1% percent by vol error
R(2,2) = (0.002/2)^2; % D sensor has 0.2% by vol error
R(3,3) = (0.2/2)^2; % internal thermistor is +/- 0.2deg
R(4,4) = (1/2)^2; % external thermistor is +/- 1deg

xs = zeros(tsteps,length(x0));
ys = zeros(tsteps-1,4);
xs(1,:) = x0;

for i = 2:tsteps
    xNoise = sqrt(Q)*randn(length(x0),1); 
    yNoise = sqrt(R)*randn(4,1);
    muC = mu_max*Ki*xs(i-1,1)/(Ki+xs(i-1,1))/(Ks+xs(i-1,1))*xs(i-1,2) ;
    xs(i,:) = f(xs(i-1,:).',S0,Yxs,m,Fs(i-1,:),Kq,V,muC,Pp,Ip(i-1,:),Rp,Np,N_power,hA,w(i-1,:),C_heat,dt) + xNoise;
    kLa = (w(i-1)^1.7)/(Fg(i-1)^0.1);
    ys(i-1,:) = g(xs(i,:).',mf_in,u(i-1,4),kLa,pr);
end

%% EKF
mu0 = x0;
cov0 = 0.001*eye(length(mu0));
mus = zeros(tsteps,length(mu0));
covs = zeros(tsteps,length(mu0),length(mu0));
mus(1,:) = mu0;
covs(1,:,:) = cov0;
for i = 2:tsteps
   [muout, covout] = Monod_EKF(mus(i-1,:).', squeeze(covs(i-1,:,:)), ys(i-1,:).', u(i-1,:), Q,R, mu_max, Yxs, m, Ki, Ks, V, S0, mf_in, Kq, Pp, Rp, Np, hA, C_heat, pr,N_power, dt);
    mus(i,:) = muout;
    covs(i,:,:) = covout;
end

%% Plotting
figure('Name','Concentrations'); hold on;
ylabel('Concentration, g/L'); xlabel('time, sec');
plot(tspan, xs(:,1),'DisplayName','Substrate Truth');
plot(tspan, xs(:,2),'DisplayName','Cell Truth');
plot(tspan, mus(:,1),'--','DisplayName','Substrate EKF');
plot(tspan, mus(:,2),'--','DisplayName','Cell EKF');
legend();

figure('Name', 'O2 Consumption'); hold on;
title('Specific O_2 Consumption Rate ~ Q_{O2}');
ylabel('Specific Consumption, mol O_2 per sec per hour'); xlabel('time, sec');
plot(tspan, xs(:,3),'DisplayName','Q_{O2} Truth');
plot(tspan, mus(:,3),'DisplayName','Q_{O2} EKF');
legend();

figure('Name', 'OR'); hold on;
title('Oxygen Transfer Rate = Oxygen Uptake Rate');
ylabel('O_R, moles O2/sec'); xlabel('time, sec');
plot(tspan, xs(:,2).*xs(:,3),'DisplayName','O_R Truth');
plot(tspan, mus(:,2).*mus(:,3),'DisplayName','O_R EKF');
legend();

figure('Name', 'Temps'); hold on;
title('System Temperatures');
ylabel('Temperature, deg F'); xlabel('time, sec');
plot(tspan,convtemp(xs(:,4),'K','F'),'DisplayName','T_C Truth');
plot(tspan,convtemp(xs(:,5),'K','F'),'DisplayName','T_\infty Truth');
plot(tspan,convtemp(mus(:,4),'K','F'),'DisplayName','T_C EKF');
plot(tspan,convtemp(mus(:,5),'K','F'),'DisplayName','T_\infty EKF');
legend();

%% Helper functions
function y = g(x,mf_in,Fg,kLa,pr)
        y = zeros(4,1);
        y(1) = mf_in - x(3)*x(2)/Fg; % Measured mass fraction of O2 leaving system
        y(2) = 1 - 1/mf_in*(1/Fg - 1/kLa/pr)*x(3)*x(2); % Measured DO2
        y(3) = x(4); % Measured Tc
        y(4) = x(5); % measured Tinf
end
function xt = f(x,S0,Yxs,m,Fs,Kq,V,muC,Pp,Ip,Rp,Np,N_power,hA,w,C_heat,dt)
   xt = zeros(size(x));

   xt(1) = x(1) + dt*(S0*Fs/V - Fs/V*x(1) - muC/Yxs - m*x(2)); % S = substrate conc
   xt(2) = x(2) + dt*(-Fs/V*x(2) + muC); % C = cell conc
   xt(3) = x(3); % QO2 = OUR coefficient
   xt(4) = x(4) + dt/C_heat*(N_power*(w^3)+Kq*x(3)*x(2)+(Pp*Ip+(Ip^2)*Rp)*Np + hA*(x(4)-x(5))); % Tc = culture temp
   xt(5) = x(5); % Tinf = environment temp
end