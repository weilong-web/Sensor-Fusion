function [xs, Ps, xf, Pf, xp, Pp] = nonLinRTSsmoother(Y, x_0, P_0, f, Q, h, R, sigmaPoints, type)
%NONLINRTSSMOOTHER Filters measurement sequence Y using a 
% non-linear Kalman filter. 
% your code here!
N = size(Y,2);
n = length(x_0);
xf = zeros(n,N);
Pf = zeros(n,n,N);
xp = zeros(n,N);
Pp = zeros(n,n,N);
xf(:,1) = x_0;
Pf(:,:,1) = P_0;
[xf, Pf, xp, Pp] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q, h, R, type);
    
% initialize outputs
xs(:,N) = xf(:,N);
Ps(:,:,N) = Pf(:,:,N);


    
% for i = 1:N
%     [xp(:,i), Pp(:,:,i)] = nonLinKFprediction(xf(:, i), Pf(:, :, i), f, T, Q, sigmaPoints, type);
%     [xf(:, i+1), Pf(:, :, i+1)] = nonLinKFupdate(xp(:,i), Pp(:,:,i), Y(:, i), S(:, i), h, R, sigmaPoints, type);
%  end
%  xf = xf(:,2:end);
%  Pf = Pf(:,:,2:end);
for k=N-1:-1:1
    [xs(:,k), Ps(:,:,k)] = nonLinRTSSupdate(xs(:,k+1), Ps(:,:,k+1), xf(:,k), Pf(:,:,k), xp(:,k+1),  Pp(:,:,k+1), f, sigmaPoints, type);
end 
%  % RTS smoother
% xs(:, N) = xf(:, N);
% Ps(:, :, N) = Pf(:, :, N);
% for i = N-1:-1:1
%     [xs(:, i), Ps(:, :, i)] = nonLinRTSSupdate(xs(:, i+1), Ps(:, :, i+1), xf(:, i), Pf(:, :, i), xp(:, i+1), Pp(:, :, i+1), f, T, sigmaPoints, type);
% end

end

function [xs, Ps] = nonLinRTSSupdate(xs_kplus1, Ps_kplus1, xf_k, Pf_k, xp_kplus1, Pp_kplus1, ...
                                     f, sigmaPoints, type)
Pk_kp1 = zeros(size(xs_kplus1,1));
if strcmp(type, 'EKF')
    [f_value,df_k] = f(xf_k);
    Pk_kp1 = Pf_k*df_k';
else
    [sp,w] = sigmaPoints(xf_k,Pf_k,type);
    for i = 1:numel(w)
        Pk_kp1 = Pk_kp1 + (sp(:,i)-xf_k)*(f(sp(:,i))-xp_kplus1).' * w(i);
    end
end
xs = xf_k + Pk_kp1 * Pp_kplus1^(-1) * (xs_kplus1 - xp_kplus1);
Ps = Pf_k - Pk_kp1 * Pp_kplus1^(-1) * (Pp_kplus1 - Ps_kplus1) * (Pk_kp1 * Pp_kplus1^(-1))';   
end
