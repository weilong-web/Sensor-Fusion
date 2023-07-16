function [xs, Ps] = nonLinRTSSupdate(xs_kplus1, ...
                                     Ps_kplus1, ...
                                     xf_k, ... 
                                     Pf_k, ...
                                     xp_kplus1, ...
                                     Pp_kplus1, ...
                                     f, ...
                                     T, ...
                                     sigmaPoints, ...
                                     type)
%NONLINRTSSUPDATE Calculates mean and covariance of smoothed state
% density, using a non-linear Gaussian model.
%
%Input:
%   xs_kplus1   Smooting estimate for state at time k+1
%   Ps_kplus1   Smoothing error covariance for state at time k+1
%   xf_k        Filter estimate for state at time k
%   Pf_k        Filter error covariance for state at time k
%   xp_kplus1   Prediction estimate for state at time k+1
%   Pp_kplus1   Prediction error covariance for state at time k+1
%   f           Motion model function handle
%   T           Sampling time
%   sigmaPoints Handle to function that generates sigma points.
%   type        String that specifies type of non-linear filter/smoother
%
%Output:
%   xs          Smoothed estimate of state at time k
%   Ps          Smoothed error convariance for state at time k

% Your code here.
Pk_kp1 = zeros(size(xs_kplus1,1));
if strcmp(type, 'EKF')
    [f_value,df_k] = f(xf_k,T);
    Pk_kp1 = Pf_k*df_k';
else
    [sp,w] = sigmaPoints(xf_k,Pf_k,type);
    for i = 1:numel(w)
        Pk_kp1 = Pk_kp1 + (sp(:,i)-xf_k)*(f(sp(:,i),T)-xp_kplus1).' * w(i);
    end
end
xs = xf_k + Pk_kp1 * Pp_kplus1^(-1) * (xs_kplus1 - xp_kplus1);
Ps = Pf_k - Pk_kp1 * Pp_kplus1^(-1) * (Pp_kplus1 - Ps_kplus1) * (Pk_kp1 * Pp_kplus1^(-1))';

    
end