function [x, P] = mu_g(x, P, yacc, Ra, g0)
% Calculate predicted accelerometer measurement using the current state estimate
hx = Qq(x)'*g0;

% Calculate the Jacobian matrix of the predicted accelerometer measurement with respect to the quaternion
[dq1, dq2, dq3, dq4] = dQqdq(x);
dhx = [dq1'*g0 dq2'*g0 dq3'*g0 dq4'*g0];
% Calculate the innovation covariance matrix
sk = dhx*P*dhx'+Ra;
% Calculate the Kalman gain
kk = P*dhx'/sk;
% Update the state estimate using the Kalman gain and accelerometer measurement
x = x+kk*(yacc-hx);
% Update the covariance matrix
P = P-kk*sk*kk';

end

