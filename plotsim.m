function plotsim(constants, tout, x, mu, cov, y, u, targ, targ_err)
    %% For plotting output of simulation
   figure('Name','Concentrations'); hold on;
ylabel('Concentration, g/L'); xlabel('time, sec');
cov_patch(1,'Substrate 95% CI');
cov_patch(2,'Cell 95% CI');
plot(tout, x(:,1),'DisplayName','Substrate Truth');
plot(tout, x(:,2),'DisplayName','Cell Truth');
plot(tout, mu(:,1),'DisplayName','Substrate EKF');
plot(tout, mu(:,2),'--','DisplayName','Cell EKF');
legend();

figure('Name', 'O2 Consumption'); hold on;
title('Specific O_2 Consumption Rate ~ Q_{O2}');
ylabel('Specific Consumption, mol O_2 per sec per hour'); xlabel('time, sec');
plot(tout, x(:,3),'DisplayName','Q_{O2} Truth');
plot(tout, mu(:,3),'DisplayName','Q_{O2} EKF');
legend();

figure('Name', 'OR'); hold on;
title('Oxygen Transfer Rate = Oxygen Uptake Rate');
ylabel('O_R, moles O2/sec'); xlabel('time, sec');
plot(tout, x(:,2).*x(:,3),'DisplayName','O_R Truth');
plot(tout, mu(:,2).*mu(:,3),'DisplayName','O_R EKF');
legend();

figure('Name', 'Temps'); hold on;
title('System Temperatures');
ylabel('Temperature, deg F'); xlabel('time, sec');
plot(tout,convtemp(x(:,4),'K','F'),'DisplayName','T_C Truth');
plot(tout,convtemp(x(:,5),'K','F'),'DisplayName','T_\infty Truth');
plot(tout,convtemp(mu(:,4),'K','F'),'DisplayName','T_C EKF');
plot(tout,convtemp(mu(:,5),'K','F'),'DisplayName','T_\infty EKF');
legend();


    function cov_patch(idx, descrip)
       myspread = [mu(:,idx)+2*sqrt(squeeze(cov(idx,idx,:))); flip(mu(:,idx)-2*sqrt(squeeze(cov(idx,idx,:))))];
       patch([tout; flip(tout)], myspread , [0.4 0.4 0.4], 'FaceAlpha',0.1,'EdgeColor','none','DisplayName',descrip);
    end
end