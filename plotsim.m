function plotsim(constants, tout, x, mu, cov, y, u, targ, targ_err)
    %% SIM OUTPUT
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
    cov_patch(3,'Q_{O2} 95% CI');
    plot(tout, x(:,3),'DisplayName','Q_{O2} Truth');
    plot(tout, mu(:,3),'DisplayName','Q_{O2} EKF');
    legend();
    
    % Calculate OR mean and cov
    mu_OR = mu(:,2).*mu(:,3);
    cov_OR = squeeze(cov(2,2,:)).*(mu(:,3).^2) + squeeze(cov(3,3,:)).*(mu(:,2).^2);
    
    figure('Name', 'OR'); hold on;
    title('Oxygen Transfer Rate = Oxygen Uptake Rate');
    ylabel('O_R, moles O2/sec'); xlabel('time, sec');
    ORspread = [mu_OR+2*sqrt(cov_OR); flip(mu_OR-2*sqrt(cov_OR))];
    patch([tout; flip(tout)], ORspread , [0.4 0.4 0.4], 'FaceAlpha',0.1,'EdgeColor','none','DisplayName','O_R 95% CI');
    plot(tout, x(:,2).*x(:,3),'DisplayName','O_R Truth');
    plot(tout, mu_OR,'DisplayName','O_R EKF');
    legend();

    figure('Name', 'Temps'); hold on;
    title('System Temperatures');
    ylabel('Temperature, deg F'); xlabel('time, sec');
    cov_patch(4,'T_C 95% CI');
    cov_patch(5,'T_\infty 95% CI');
    plot(tout,convtemp(x(:,4),'K','F'),'DisplayName','T_C Truth');
    plot(tout,convtemp(x(:,5),'K','F'),'DisplayName','T_\infty Truth');
    plot(tout,convtemp(mu(:,4),'K','F'),'DisplayName','T_C EKF');
    plot(tout,convtemp(mu(:,5),'K','F'),'DisplayName','T_\infty EKF');
    legend();

    %% SIM CONTROL
    figure('Name', 'Substrate Feed'); hold on;
    title('Substrate Feed Input');
    ylabel('\nu_S, L/s'); xlabel('time, sec');
    plot(tout, u(:,2));
    
    figure('Name', 'Peltier'); hold on;
    title('Peltier Current Input');
    ylabel('I_p, Amps'); xlabel('time, sec');
    plot(tout, u(:,3));
    
    %% SIM TARGETS
    
    function cov_patch(idx, descrip)
       myspread = [mu(:,idx)+2*sqrt(squeeze(cov(idx,idx,:))); flip(mu(:,idx)-2*sqrt(squeeze(cov(idx,idx,:))))];
       if idx > 3
          myspread = convtemp(myspread,'K','F'); 
       end
       patch([tout; flip(tout)], myspread , [0.4 0.4 0.4], 'FaceAlpha',0.1,'EdgeColor','none','DisplayName',descrip);
    end
end