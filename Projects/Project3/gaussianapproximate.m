function [mu_y, Sigma_y, y_s, x_s] = gaussianapproximate(mu_x, Sigma_x, f, N)

    if nargin < 4
        N = 5000;
    end

    % sample in the original gaussian distribution
    x_s = mvnrnd(mu_x, Sigma_x, N)';
    % apply general non-linear transformation function to samples
    y_s = f(x_s);
    % calculate mean of the transformed samples
    mu_y = mean(y_s,2);
    % calculate estimated unbiased covariance of the transformed samples
    Sigma_y = 1/(N-1) * (y_s - mu_y) * (y_s - mu_y)';
end