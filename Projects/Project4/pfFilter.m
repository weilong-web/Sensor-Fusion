function [xfp, Pfp, Xp, Wp] = pfFilter(x_0, P_0, Y, proc_f, proc_Q, meas_h, meas_R, ...
                             N, bResample, plotFunc)
%PFFILTER Filters measurements Y using the SIS or SIR algorithms and a
% state-space model.
%
% Input:
%   x_0         [n x 1] Prior mean
%   P_0         [n x n] Prior covariance
%   Y           [m x K] Measurement sequence to be filtered
%   proc_f      Handle for process function f(x_k-1)
%   proc_Q      [n x n] process noise covariance
%   meas_h      Handle for measurement model function h(x_k)
%   meas_R      [m x m] measurement noise covariance
%   N           Number of particles
%   bResample   boolean false - no resampling, true - resampling
%   plotFunc    Handle for plot function that is called when a filter
%               recursion has finished.
% Output:
%   xfp         [n x K] Posterior means of particle filter
%   Pfp         [n x n x K] Posterior error covariances of particle filter
%   Xp          [n x N x K] Non-resampled Particles for posterior state distribution in times 1:K
%   Wp          [N x K] Non-resampled weights for posterior state x in times 1:K

% Your code here, please. 
% If you want to be a bit fancy, then only store and output the particles if the function
% is called with more than 2 output arguments.
    % Extract dimensions
    n = size(x_0, 1);  % State dimension
    K = size(Y, 2);    % Number of time steps

    % Initialize variables
    xfp = zeros(n, K);
    Pfp = zeros(n, n, K);
    Xp = zeros(n, N, K);
    Wp = zeros(N, K);

    % Initialize prior
    Xk = mvnrnd(x_0, P_0, N)';
    Wk = ones(1, N) / N;

    % Filter loop
    for k = 1:K
        yk = Y(:, k);

        % Particle filter step
        [Xk, Wk] = pfFilterStep(Xk, Wk, yk, proc_f, proc_Q, meas_h, meas_R);

        % Store filtered particles and weights
        Xp(:, :, k) = Xk;
        Wp(:, k) = Wk;

        % Compute filtered estimate and covariance
        xfp(:, k) = Xk * Wk';
        Pfp(:, :, k) = (Xk - xfp(:, k)) * diag(Wk) * (Xk - xfp(:, k))';

        % Resampling
        if bResample
            [Xk, Wk, ~] = resampl(Xk, Wk);
        end
   
    end
end

function [Xk, Wk, j] = resampl(X, W)

    N=size(X,2);
    
    % Generates the segmented numberline from 0 to 1
    segment = [0 cumsum(W)/sum(W)];
    
    % draw samples from uniform distribution on [0,1]
    samples = rand([1 N]);
    
    j=zeros(1,N);
    for i=1:N
        j(i) = find(samples(i) >= segment,1,'last');
    end
    
    Wk = 1/N*ones(1,N);
    Xk = X(:,j);
end

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