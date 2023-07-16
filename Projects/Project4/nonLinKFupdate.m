function [x, P] = nonLinKFupdate(x, P, y, h, R, type)
%NONLINKFUPDATE calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   y           [m x 1] measurement vector
%   h           Measurement model function handle
%               [hx,Hx]=h(x) 
%               Takes as input x (state), 
%               Returns hx and Hx, measurement model and Jacobian evaluated at x
%               Function must include all model parameters for the particular model, 
%               such as sensor position for some models.
%   R           [m x m] Measurement noise covariance
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] updated state mean
%   P           [n x n] updated state covariance
%
    n = size(x,1);
    switch type
        case 'EKF'
            
        [hx, dhx] = h(x);
        %inovation covariance and kalman gain
        S = dhx * P * dhx' + R;
        K = P * dhx' / S;
        % update
        x = x + K * ( y - hx );
        P = P - K * S * K';   
            
        case 'UKF'
    
        %sigma points
        [SP,W] = sigmaPoints(x, P, type);
        yhat = 0*y;
        for i=1:numel(W)
            yhat = yhat + h(SP(:,i)) * W(i);
        end
        Pxy = 0*x*y';
        for i=1:numel(W)
            Pxy = Pxy + (SP(:,i)-x)*(h(SP(:,i))-yhat)' * W(i);
        end
        S = R;
        for i=1:numel(W)
            S = S + (h(SP(:,i))-yhat)*(h(SP(:,i))-yhat)' * W(i);
        end

        %mean and covariance
        x = x + Pxy / S * ( y - yhat );
        P = P - Pxy / S * Pxy';

        % Make sure the covariance matrix is semi-definite
        if min(eig(P))<=0
            [v,e] = eig(P, 'vector');
            e(e<0) = 1e-4;
            P = v*diag(e)/v;
        end
            
        case 'CKF'
    
        % compute sigma points
        [SP,W] = sigmaPoints(x, P, type);

        yhat = 0*y;
        for i=1:numel(W)
            yhat = yhat + h(SP(:,i)) * W(i);
        end
        Pxy = 0*x*y';
        for i=1:numel(W)
            Pxy = Pxy + (SP(:,i)-x)*(h(SP(:,i))-yhat)' * W(i);
        end
        S = R;
        for i=1:numel(W)
            S = S + (h(SP(:,i))-yhat)*(h(SP(:,i))-yhat)' * W(i);
        end

        %updated mean and covariance
        x = x + Pxy / S * ( y - yhat );
        P = P - Pxy / S * Pxy';

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

