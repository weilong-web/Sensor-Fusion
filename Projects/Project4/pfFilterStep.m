function [X_k, W_k] = pfFilterStep(X_kmin1, W_kmin1, yk, proc_f, proc_Q, meas_h, meas_R)
    N = size(X_kmin1, 2);  % Number of particles
    n = size(X_kmin1, 1);  % State dimension
    m = size(yk, 1);       % Measurement dimension
    X_k = mvnrnd(proc_f(X_kmin1)' ,proc_Q )';  
    
    % calculate p(y_k|x(i)_k) 
    Wy = mvnpdf(yk',meas_h(X_k)',meas_R)';
    
    % compute weights
    W_k = W_kmin1 .* Wy;
    W_k = W_k / sum(W_k);
end