%% pH Test
clearvars; close all; clc;

%% Parameters
dt = 1;
tend = 3600;
tspan = 0:dt:tend;
R = (0.1/2)^2 ; % 2sigma confidence interval of pH is +/- 0.1
Q = 1E-14*dt^2*eye(3); % process noise
rN0 = 0; % NH3 consumption rate of bacteria in moles/liter/sec
pH0 = 7;
x0 = [10^(-pH0); 0; rN0]; % [H30+, NH4+, rNH3] in [moles/liter, moles/liter, moles/sec]

%% Control Inputs
CN = 0.2; % concentration percentage of NH4OH moles/liter in feed line (assumes 15% ammonia in air)
FN = 0*[0.01*ones(1,round(length(tspan)/3+1,0)), 0*ones(1,round(length(tspan)/3+1,0)), 0.02*ones(1,round(length(tspan)/3+1,0))].'; 
% volumetric flow rate of pH feed, liters per sec
V = 5; % volume of bioreactor, liters
Ka = 5.6e-10; % acid disassociation constant for NH4+

%% Simulation
x = zeros(length(tspan),length(x0)); % ground truth
x(1,:) = x0;
y = zeros(length(tspan)-1,1); % one measurement, pH
for i = 2:length(tspan)
    Vnoise = sqrt(R)*randn(1);
    Wnoise = sqrt(Q)*randn(3,1);
    
    NH3 = Ka*x(i-1,2)/x(i-1,1) ;
    deltaNH3 = -dt*x(i-1,3) + dt*CN*FN(i-1)/V ;
    mysign = sign(deltaNH3);
    if mysign == 0
       mysign = -1; 
    end
    NH3 = NH3+ deltaNH3 ;
    b = x(i-1,1) + NH3 + Ka;
    c = x(i-1,1) * NH3 - x(i-1,2)* Ka;
    X = 0.5*(b + mysign*sqrt(b^2 - 4*c));
    x(i,1) = x(i-1,1) - X + Wnoise(1);
    x(i,2) = x(i-1,2) + X + Wnoise(2);
    x(i,3) = x(i-1,3) + Wnoise(3);
    y(i-1) = -log10(x(i,1)) + Vnoise;
end
pHtrue = -log10(x(:,1));
%% EKF
C =@(x,FN) [-log10(exp(1))/x(1) 0 0] ;
A =@(x,FN) myA(x,FN,dt,CN,V,Ka);
f =@(x,FN) myf(x,FN,dt,CN,V,Ka);
g =@(x,FN) -log10(x(1));
mu0 = x0;
cov0 = 0.1*eye(3);

[mu, cov] = EKF(mu0, cov0, tspan, y, f, g, FN, A, C, Q, R);
pHEKF = -log10(mu(:,1)); % calculate pH from [H30+]
%% Plot
figure('Name','pH'); hold on;
xlabel('time, sec'); ylabel('pH');
plot(tspan, pHtrue,'DisplayName','True');
plot(tspan, pHEKF,'DisplayName','EKF Estimate');
scatter(tspan(2:end), y, 'DisplayName' ,'Measured');
legend;

%% Functions
function At = myA(x,FN,dt,CN,V,Ka)
    NH3 = Ka*x(2)/x(1) ;
    deltaNH3 = - dt*x(3) + dt*CN*FN/V ;
    mysign = sign(deltaNH3);
    if mysign == 0
       mysign = -1; 
    end
    NH3 = NH3+ deltaNH3 ;
    At = eye(3);
    b = x(1) + NH3 + Ka;
    c = x(1) * NH3 - x(2)* Ka;
    disc = sqrt(b^2 - 4*c);
    At(2,1) = 0.5*(1-Ka*x(2)/(x(3)^2)+mysign*(x(1)-(Ka^2)*(x(2)^2)/(x(1)^3)+deltaNH3+Ka)/disc );
    At(1,1) = 1-At(2,1);
    At(1,2) = -0.5*(Ka/x(1) + mysign*((deltaNH3+Ka)*Ka/x(1) + (Ka^2)*x(2)/x(1) + Ka)/disc);
    At(2,2) = 1-At(1,2);
    At(1,3) = -0.5*(-dt+mysign*((dt^2)*(x(3)-CN*FN/V)+2*dt*x(1))/disc) ;
    At(2,3) = -At(1,3) ;
end
function ft = myf(x,FN,dt,CN,V,Ka)
    ft = zeros(size(x));
    NH3 = Ka*x(2)/x(1) ;
    deltaNH3 = - dt*x(3) + dt*CN*FN/V ;
    mysign = sign(deltaNH3);
    if mysign == 0
       mysign = -1; 
    end
    NH3 = NH3+ deltaNH3 ;
    b = x(1) + NH3 + Ka;
    c = x(1) * NH3 - x(2)* Ka;
    X = 0.5*(b + mysign*sqrt(b^2 - 4*c));
    ft(1) = x(1) - X;
    ft(2) = x(2) + X;
    ft(3) = x(3);
end
