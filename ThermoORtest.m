%% With dOR
%Use SI measurements - J, s, m, L, 

clear
close all

%% Parameters to edit

dt = 5;
t = 0:dt:(24*60*60);
n = length(t);

%Constants
Cc = 836800; 
N = 2*10^-10;
Kb = 470000; % bacteria heat coefficient, 
Pp = 20; % peltier coefficient, W/A
Rp = 2; % peltier resistance, ohms
hA = 600*.0001;

%OR model
OR1 = 0;
OR2 = .0001;
ttrans = 12*60*60;
width = 4*60*60;

%Tinf model
T = 24*60*60;
base = 21;
delta = 2;

%Current model
Tcurr = 30*60;
Ip = 2.5*cos(t*2*pi/Tcurr)-2.5;

%Control inputs
rpm = [500*ones(1,round(n/3+1,0)), 3500*ones(1,round(n/3+1,0)), 750*ones(1,round(n/3+1,0))];
gas_flow = (50/1000)*(1 + t./(8*3600))/60; %sccm to slps

%% Model

%states = [Tc Tinf OR dOR]
%measurements = [Tc Tinf mfO2out %DO]

R = diag([.2^2, .5^2, .002^2, .001^2]); %sensor noise
q = (.000001)^2;  %process noise
Q = diag([(dt*.002)^2, (dt*.002)^2, q, 4*q/(dt^2)]);

x(:,1) = [23;21;(OR2 - OR1)/(1 + exp(-(t(1)-ttrans)/width)) + OR1; 0];

for i = 1:length(t)

    V = chol(R)*randn(4,1);
    W = chol(Q)*randn(4,1);
    
    kLa(i) = .2*(rpm(i)^1.7)*(gas_flow(i)^-0.1)/(1000*60*500); %
    
    x(1,i+1) = x(1,i) + dt*(N*(rpm(i)^3) + Kb*x(3,i) + Pp*Ip(i) + (Ip(i)^2)*Rp - hA*(x(1,i)-x(2,i)))/Cc + W(1);
    x(2,i+1) = delta*sin(t(i)*2*pi/T) + base + W(2);
    x(3,i+1) = (OR2 - OR1)/(1 + exp(-(t(i)-ttrans)/width)) + OR1 + W(3);
    x(4,i+1) = (OR2 - OR1)*exp(-(t(i)-ttrans)/width)/(width*(exp(-(t(i)-ttrans)/width)+1)^2);% + W(2);
    y(:,i) = [x(1,i); x(2,i);.21 - x(3,i)*(1/gas_flow(i)); ...
        1 - x(3,i)*(1/gas_flow(i) + 1/kLa(i))/.21] + V; 
    z(1,i+1) = N*(rpm(i)^3);
    z(2,i+1) = Kb*x(3,i);
    z(3,i+1) = Pp*Ip(i) + (Ip(i)^2)*Rp ;
    z(4,i+1) = - hA*(x(1,i)-x(2,i));
    z(1:4,i+1) = dt*z(1:4,i+1)/Cc;
    
end

x = x(:,1:(end-1));
z = z(:,1:(end-1));


figure()
plot(t,z(1,:),t,z(2,:),t,z(3,:),t,z(4,:))
xlabel("Time (s)")
ylabel("Value")
title("OR")

figure()
plot(t,x(3,:))
xlabel("Time (s)")
ylabel("Value")
title("OR")

figure()
plot(t,x(4,:))
xlabel("Time (s)")
ylabel("dOR")
title("Ideal dOR")

figure()
plot(t,y(3,:))
xlabel("Time (s)")
ylabel("Value")
title("mO2")

figure()
plot(t,y(4,:))
xlabel("Time (s)")
ylabel("Value")
title("%O2")

%% EKF

clear mu sig

mu(:,1) = [23; 21; 0; 0];
sig(:,:,1) = diag([0.1,0.25,0.01,0.000001]);
time_ave = 0;
A = eye(4);
A(1,:) = [1 - dt*hA/Cc, dt*hA/Cc, dt*Kb/Cc,0];
A(3,4) = dt;
C = eye(4);

for i = 1:(length(t)-1)
    
    mu_pre = A*mu(:,i);
    sig_pre = A*sig(:,:,i)*(A') + Q;

    C(3:4,3:4) = [-1/gas_flow(i), 0; -(1/gas_flow(i) + 1/kLa(i))/.21, 0];

    y_pre = [mu_pre(1); mu_pre(2);.21 - mu_pre(3)*(1/gas_flow(i)); 1 - mu_pre(3)*(1/gas_flow(i) + 1/kLa(i))/.21];
    K = sig_pre*(C')*inv(C*sig_pre*(C') + R);
    mu(:,i+1) = mu_pre + K*(y(:,i+1) - y_pre);
    sig(:,:,i+1) = sig_pre - K*C*sig_pre;
    
    
end 

figure()
plot(t,x(1,:),t,mu(1,:))
xlabel("Time (s) ")
ylabel("Tc")
title("EKF - T_C Output")
legend(["True value", "Estimated value"],"Location","best")

figure()
plot(t,x(2,:),t,mu(2,:))
xlabel("Time (s) ")
ylabel("Tinf")
title("EKF - T_\infty Output")
legend(["True value", "Estimated value"],"Location","best")

figure()
plot(t,x(3,:),t,mu(3,:))
xlabel("Time (s) ")
ylabel("O_R")
title("EKF - O_R Output")
legend(["True value", "Estimated value"],"Location","best")

figure()
plot(t,mu(4,:),t,x(4,:))
xlabel("Time (s) ")
ylabel("dOR")
title("dOR")
legend(["Estimated dOR","Ideal dOR"],"Location","best")

