function dydt = Model2(t, y, theta)
% Function that describes the ODE-system for the model2. 
% Args:
%   t, the time 
%   y, the different states with:
%       y(1), intralcellular glucose 
%       y(2), Snf1p 
%       y(3), Suc2 (the outsignal)
%       y(4), X the cellular response to starvation 
%   theta, paramater vector where theta(0) -- theta(10) corresponds to the
%   rates k1, ..., k10. theta(11) corresponds to the inital glucose level.
% Returns:
%   The deriviatives for the system at time t which is used by the
%   ODE-solvers in Matlab. 

% Rename the parameters according to rates 
k1 = theta(1);
k2 = theta(2);
k3 = theta(3);
k4 = theta(4);
k5 = theta(5);
k6 = theta(6);
k7 = theta(7);
k8 = theta(8);
k9 = theta(9);
k10 = theta(10);


% Fix the rate in 
if t < 0.01
    rate_in = k1 * 1;
else 
    rate_in = k1 * 1 / 40;
end 

% The ode-functions 
dy1dt = rate_in - k2 * y(1) + k3 * y(4);
dy2dt = k4 / (k8 + y(1)^2) - k5 * y(2);
dy3dt = k6 * y(2)^2 - k7 * y(3)^2;
dy4dt = (k10 * y(2) + k7 * y(3)^2) / (y(1) * y(3) + k8) - k9 * y(4);

% Note that y3 is the out signal
dydt = [dy1dt, dy2dt, dy3dt dy4dt]';

end 