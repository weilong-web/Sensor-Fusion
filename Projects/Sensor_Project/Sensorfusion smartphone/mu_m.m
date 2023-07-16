function [x, P] = mu_m(x, P, mag, m0,Rm)

hx = Qq(x)'*m0;
[dq1, dq2, dq3, dq4] = dQqdq(x);
dhx = [dq1'*m0 dq2'*m0 dq3'*m0 dq4'*m0];

sk = dhx*P*dhx'+Rm;
kk = P*dhx'/sk;

x = x+kk*(mag-hx);
P = P-kk*sk*kk';

end