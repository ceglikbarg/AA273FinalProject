function [mu_out, cov_out] = Monod_EKF(mu_in, cov_in, y, u, Q,R, mu_max, Yxs, m, Ki, Ks, V, S0, mf_in, Kq, Pp, Rp, Np, hA, C_heat, pr,N_power, dt)
    % Collect inputs
%     S = mu_in(1); % cell concentration, g/l
%     C = mu_in(2); % substrate concentration, g/l
%     QO2 = mu_in(3); % OTR = OUR, moles/sec
%     Tc = mu_in(4); % Temperature of culture
%     Tinf = mu_in(5); % Temperature of environment
    mygrowthrate = mu(mu_in(1),mu_in(2));
    
    % Collect controls
    w = u(1); % agitator speed, rpm
    Fs = u(2); % substrate feed rate, L/s
    Ip = u(3); % Peltier current, amps
    Fg = u(4); % air feed rate, SLPS
    kLa = w^(1.7)/(Fg^0.1); % kLa coeff
    
    % Collect measurements
%     mf_out = y(1);
%     DO = y(2);
%     Tc_meas = y(3);
%     Tinf_meas = y(4);
    
    % EKF Update
    mu_p = f(mu_in);
    At = A(mu_in);
    cov_p = At*cov_in*(At.') + Q;
    
    % EKF Prediction
    Ct = C(mu_p);
    ypred = g(mu_p);
    Kt = cov_p*(Ct.')/(Ct*cov_p*(Ct.')+R); % Kalman gain
    mu_out = mu_p + Kt*(y - ypred);
    cov_out = (eye(length(mu_in)) - Kt*Ct)*cov_p;
    
    function xt = f(x)
       xt = zeros(size(x));
       
       xt(1) = x(1) + dt*(S0*Fs/V - Fs/V*x(1) - mygrowthrate/Yxs - m*x(2)); % S = substrate conc
       xt(2) = x(2) + dt*(-Fs/V*x(2) + mygrowthrate); % C = cell conc
       xt(3) = x(3); % QO2 = OUR coefficient
       xt(4) = x(4) + dt/C_heat*(N_power*(w^3)+Kq*x(3)*x(2)+(Pp*Ip+(Ip^2)*Rp)*Np + hA*(x(4)-x(5))); % Tc = culture temp
       xt(5) = x(5); % Tinf = environment temp
    end
    function At = A(x)
        At = zeros(length(x));
        mymu = mu(x(1),x(2));
        dmu = mu_max*x(2)*Ki*(Ki*Ks-(x(1)^2))/((Ki+x(1))^2)/((Ks+x(1))^2);
        At(1,1) = 1 - dt*Fs/V - dt*dmu/Yxs;
        At(1,2) = -dt*mymu - m;
        At(2,1) = dt*dmu;
        At(2,2) = 1 - dt*Fs/V + dt*mymu ;
        At(3,3) = 1;
        At(4,2) = dt/C_heat*Kq*x(3);
        At(4,3) = dt/C_heat*Kq*x(2);
        At(4,4) = 1 + dt/C_heat*hA;
        At(4,5) = -dt/C_heat*hA;
        At(5,5) = 1;
    end
    function ypred = g(x)
        ypred = zeros(4,1);
        ypred(1) = mf_in - x(3)*x(2)/Fg; % Measured mass fraction of O2 leaving system
        ypred(2) = 1 - 1/mf_in*(1/Fg - 1/kLa/pr)*x(3)*x(2); % Measured DO2
        ypred(3) = x(4); % Measured Tc
        ypred(4) = x(5); % measured Tinf
    end
    function Ct = C(x)
       Ct = zeros(4,5);
       Ct(1,2) = -x(3)/Fg;
       Ct(1,3) = -x(2)/Fg;
       Ct(2,2) = -1/mf_in * (1/Fg - 1/pr/kLa)*x(3);
       Ct(2,3) = -1/mf_in * (1/Fg - 1/pr/kLa)*x(2);
       Ct(3,4) = 1;
       Ct(4,5) = 1;
    end
    function mut = mu(S,C)
        mut = mu_max*Ki*S/(Ki+S)/(Ks+S)*C ;
    end
end

