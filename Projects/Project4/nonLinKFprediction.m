function [x, P] = nonLinKFprediction(x, P, f, Q, type)
%NONLINKFPREDICTION calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   f           Motion model function handle
%               [fx,Fx]=f(x) 
%               Takes as input x (state), 
%               Returns fx and Fx, motion model and Jacobian evaluated at x
%               All other model parameters, such as sample time T,
%               must be included in the function
%   Q           [n x n] Process noise covariance
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] predicted state mean
%   P           [n x n] predicted state covariance
%
n = size(x,1);
    switch type
        case 'EKF'
            
        [fx, dfx] = f(x);
        %Linearize
        x = fx;
        P = dfx * P * dfx' + Q;
            
        case 'UKF'
   
        [SP,W] = sigmaPoints(x, P, type);
        %mean
        x = zeros(n,1);
        for i=1:numel(W)
            x = x + f(SP(:,i)) * W(i);
        end
        %covariance
        P = Q;
        for i=1:numel(W)
            P = P + (f(SP(:,i))-x)*(f(SP(:,i))-x).' * W(i);
        end
        % Make sure the covariance matrix is semi-definite
        if min(eig(P))<=0
            [v,e] = eig(P, 'vector');
            e(e<0) = 1e-4;
            P = v*diag(e)/v;
        end
                
        case 'CKF'
            
        %sigma points
        [SP,W] = sigmaPoints(x, P, type);
        %mean
        x = zeros(n,1);
        for i=1:numel(W)
            x = x + f(SP(:,i)) * W(i);
        end
        %covariance
        P = Q;
        for i=1:numel(W)
            P = P + (f(SP(:,i))-x)*(f(SP(:,i))-x).' * W(i);
        end      
        % Make sure the covariance matrix is semi-definite
        if min(eig(P))<=0
            [v,e] = eig(P, 'vector');
            e(e<0) = 1e-4;
            P = v*diag(e)/v;
        end    
        otherwise
            error('Incorrect type of non-linear Kalman filter')
    end

end