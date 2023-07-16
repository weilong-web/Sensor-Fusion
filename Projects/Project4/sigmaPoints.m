function [SP,W] = sigmaPoints(x, P, type)
% SIGMAPOINTS computes sigma points, either using unscented transform or
% using cubature.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%
%Output:
%   SP          [n x 2n+1] UKF, [n x 2n] CKF. Matrix with sigma points
%   W           [1 x 2n+1] UKF, [1 x 2n] UKF. Vector with sigma point weights 
%

    switch type        
        case 'UKF'
            n = size(x,1);
            SP = zeros(n, 2*n+1);
            Psqrt = sqrtm(P);
            W0 = 1 - n/3;
            SP(:,1) = x;
            for i=1:n
                SP(:, i+1  ) = x + sqrt( n / (1 - W0) ) * Psqrt(:,i);
                SP(:, i+n+1) = x - sqrt( n / (1 - W0) ) * Psqrt(:,i);
            end
            W  = [W0, (1-W0)/(2*n)*ones(1, 2*n)];
                
        case 'CKF'
            n = size(x,1);           
            SP = zeros(n, 2*n);
            Psqrt = sqrtm(P);
            for i=1:n
                SP(:, i  ) = x + sqrt(n) * Psqrt(:,i);
                SP(:, i+n) = x - sqrt(n) * Psqrt(:,i);
            end
            W  = 1/(2*n) * ones(1, 2*n);
        otherwise
            error('Incorrect type of sigma point')
    end

end