%% With dOR

clear
close all

%% Parameters to edit

dt = 5;
t = 0:dt:3600;
N = length(t);

%OR model
OR1 = 0;
OR2 = 1.2;
ttrans = 1800;
width = 200;

%Control inputs
rpm = [500*ones(1,round(N/3+1,0)), 3500*ones(1,round(N/3+1,0)), 750*ones(1,round(N/3+1,0))];
gas_flow = (50/1000)*(1 + t./1800); %sccm to slpm

%% Model

R = [.002, 0; 0, .001]; %sensor noise
q = .0001;  %process noise
Q = [q, 0; 0, 2*q/dt];

x(:,1) = [(OR2 - OR1)/(1 + exp(-(t(1)-ttrans)/width)) + OR1; 0];

for i = 1:length(t)

    V = chol(R)*randn(2,1);
    W = chol(Q)*randn(2,1);
    
    kLa(i) = (rpm(i)^1.7)*(gas_flow(i)^-0.1);
      
    x(1,i+1) = (OR2 - OR1)/(1 + exp(-(t(i)-ttrans)/width)) + OR1 + W(1);
    x(2,i+1) = (OR2 - OR1)*exp(-(t(i)-ttrans)/width)/(width*(exp(-(t(i)-ttrans)/width)+1)^2);% + W(2);
    x(3,i+1) = (x(1,i+1) - x(1,i))/dt;
    y(:,i) = [.21 - x(1,i+1)*(1/(60*gas_flow(i))); ...
        1 - x(1,i+1)*(1/(60*gas_flow(i)) + 1/(kLa(i)*60))/.21] + V; 
%     if y(1,i) < 0
%         y(1,i) = 0;
%     end
%     if y(2,i) < 0
%         y(2,i) = 0;
%     end
end

x = x(:,1:(end-1));

figure()
plot(t,x(1,:))
xlabel("Time (s)")
ylabel("Value")
title("OR")

figure()
plot(t,x(2,:),t,x(3,:))
xlabel("Time (s)")
ylabel("dOR")
title("dOR")
legend(["Ideal dOR", "True dOR"],"Location","best")

figure()
plot(t,y(1,:))
xlabel("Time (s)")
ylabel("Value")
title("mO2")

figure()
plot(t,y(2,:))
xlabel("Time (s)")
ylabel("Value")
title("%O2")

clear mu sig

mu(:,1) = [0; 0];
sig(:,:,1) = 0.01*eye(2);
time_ave = 0;
A = [1, dt; 0, 1];

for i = 1:(length(t)-1)
 
    mu_pre = A*mu(:,i);
    sig_pre = A*sig(:,:,i)*(A') + Q;

    C = [-1/(60*gas_flow(i)), 0; -(1/(60*gas_flow(i)) + 1/(kLa(i)*60))/.21, 0];

    y_pre = [.21 - mu_pre(1)*(1/(60*gas_flow(i))); 1 - mu_pre(1)*(1/(60*gas_flow(i)) + 1/(kLa(i)*60))/.21];
    K = sig_pre*(C')*inv(C*sig_pre*(C') + R);
    mu(:,i+1) = mu_pre + K*(y(:,i+1) - y_pre);
    sig(:,:,i+1) = sig_pre - K*C*sig_pre;
    
    
end 

figure()
plot(t,x(1,:),t,mu(1,:))
xlabel("Time (s) ")
ylabel("OR")
title("EKF")
legend(["True value", "Estimated value"],"Location","best")

figure()
plot(t,x(2,:),t,x(3,:),t,mu(2,:))
xlabel("Time (s) ")
ylabel("dOR")
title("dOR")
legend(["Ideal dOR", "True dOR", "Estimated dOR"],"Location","best")

%% Simplified
% 
% clear
% close all
% 
% dt = 1;
% t = 0:dt:3600;
% 
% %OR model
% OR1 = 0;
% OR2 = 1.2;
% ttrans = 1800;
% width = 300;
% 
% %Control inputs
% rpm = [750*ones(1,1500), 1500*ones(1,1000), 750*ones(1,1500)];
% gas_flow = (50/1000)*(1 + t/1800); %sccm to slpm
% 
% R = [.002, 0; 0, .001];
% q = .0001;
% Q = q;
% 
% x(:,1) = [(OR2 - OR1)/(1 + exp(-(t(1)-ttrans)/width)) + OR1; 0];
% 
% for i = 1:length(t)
% 
%     V = chol(R)*randn(2,1);
%     W = chol(Q)*randn(2,1);
%     
%     kLa(i) = (rpm(i)^1.7)*(gas_flow(i)^-0.1);
%       
%     x(1,i+1) = (OR2 - OR1)/(1 + exp(-(t(i)-ttrans)/width)) + OR1 + W(1);
%     y(:,i) = [.21 - x(1,i+1)*(1/(60*gas_flow(i))); ...
%         1 - x(1,i+1)*(1/(60*gas_flow(i)) + 1/(kLa(i)*60))/.21] + V; 
% %     if y(1,i) < 0
% %         y(1,i) = 0;
% %     end
% %     if y(2,i) < 0
% %         y(2,i) = 0;
% %     end
% end
% 
% x = x(:,1:(end-1));
% 
% figure()
% plot(t,x(1,:))
% xlabel("Time (s)")
% ylabel("Value")
% %title("EKF")
% 
% figure()
% plot(t,y(1,:))
% xlabel("Time (s)")
% ylabel("Value")
% 
% figure()
% plot(t,y(2,:))
% xlabel("Time (s)")
% ylabel("Value")
% 
% clear mu sig
% 
% mu(:,1) = [0];
% sig(:,:,1) = 0.01;
% time_ave = 0;
% 
% for i = 1:(length(t)-1)
%  
%     A = 1;
%     mu_pre = mu(:,i);
%     sig_pre = A*sig(:,:,i)*(A') + Q;
% 
%     C = [-1/(60*gas_flow(i)); -(1/(60*gas_flow(i)) + 1/(kLa(i)*60))/.21];
% 
%     y_pre = [.21 - mu_pre(1)*(1/(60*gas_flow(i))); 1 - mu_pre(1)*(1/(60*gas_flow(i)) + 1/(kLa(i)*60))/.21];
%     K = sig_pre*(C')*inv(C*sig_pre*(C') + R);
%     mu(:,i+1) = mu_pre + K*(y(:,i+1) - y_pre);
%     sig(:,:,i+1) = sig_pre - K*C*sig_pre;
%     
%     
% end 
% 
% figure()
% plot(t,x(1,:),t,mu(1,:))
% xlabel("Time (s) ")
% ylabel("OR")
% title("EKF")
% legend(["True value", "Estimated value"],"Location","best")
% 
