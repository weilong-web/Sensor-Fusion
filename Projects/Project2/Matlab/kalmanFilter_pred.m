function [x, P, x_pred, P_pred] = kalmanFilter(Y, x_0, P_0, A, Q, H, R)
%KALMANFILTER Filters measurements sequence Y using a Kalman filter.
%
%Input:
% Y [m x N] Measurement sequence
% x_0 [n x 1] Prior mean
% P_0 [n x n] Prior covariance
% A [n x n] State transition matrix
% Q [n x n] Process noise covariance
% H [m x n] Measurement model matrix
% R [m x m] Measurement noise covariance
%
%Output:
% x [n x N] Estimated state vector sequence
% P [n x n x N] Filter error convariance
% x_pred [n x N] Predicted state vector sequence
% P_pred [n x n x N] Predicted error covariance sequence

%% Parameters
N = size(Y,2);

n = length(x_0);
m = size(Y,1);

%% Data allocation
x = zeros(n,N);
P = zeros(n,n,N);
x_pred = zeros(n,N);
P_pred = zeros(n,n,N);
x(:,1) = x_0;
P(:,:,1) = P_0;

for k = 1:N
% Predict step
    if k > 1
        x_pred(:,k) = A * x(:,k-1);
        P_pred(:,:,k) = A * P(:,:,k-1) * A' + Q;
    else
        x_pred(:,k) = A * x_0;
        P_pred(:,:,k) = A * P_0 * A' + Q;
    end
% Update step
    K = P_pred(:,:,k) * H' * inv(H * P_pred(:,:,k) * H' + R);
    x(:,k) = x_pred(:,k) + K * (Y(:,k) - H * x_pred(:,k));
    P(:,:,k) = (eye(n) - K * H) * P_pred(:,:,k);
end
end

