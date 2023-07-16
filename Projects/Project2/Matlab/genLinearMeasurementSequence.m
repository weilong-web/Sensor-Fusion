function Y = genLinearMeasurementSequence(X, H, R)
%GENLINEARMEASUREMENTSEQUENCE generates a sequence of observations of the state 
% sequence X using a linear measurement model. Measurement noise is assumed to be 
% zero mean and Gaussian.
%
%Input:
%   X           [n x N+1] State vector sequence. The k:th state vector is X(:,k+1)
%   H           [m x n] Measurement matrix
%   R           [m x m] Measurement noise covariance
%
%Output:
%   Y           [m x N] Measurement sequence
%

% Compute number of measurements and number of time steps
m = size(H, 1);
N = size(X, 2) - 1;

% Initialize measurement sequence
Y = zeros(m, N);

% Generate measurements using measurement model
for k = 1:N
    Y(:,k) = H*X(:,k+1) + mvnrnd(zeros(m,1), R)';
end
end