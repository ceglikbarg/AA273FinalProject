function [mu, cov] = EKF(mu0, cov0, tspan, y, f, g, u, A, C, Q, R)
    mu = zeros(length(tspan),length(mu0));
    cov = zeros(length(tspan), size(cov0,1), size(cov0,2));
    mu(1,:) = mu0;
    cov(1,:,:) = cov0;
    for i = 2:length(tspan)
        % Predict
        mu_p = f(mu(i-1,:).',u(i-1,:)) ;
        At = A(mu(i-1,:).',u(i-1,:));
        cov_p = At*(squeeze(cov(i-1,:,:)))*(At.') + Q;
        
        % Update
        Ct = C(mu_p,u(tspan(i)));
        K = cov_p*(Ct.')/(Ct*cov_p*(Ct.')+R); % Kalman gain
        mu(i,:) = mu_p + K*(y(i-1,:).' - g(mu_p,u(i,:)));
        cov(i,:,:) = cov_p - K*Ct*cov_p;
    end
end