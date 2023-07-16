function [x, P] = tu_qw(x, P, omega, T, Rw)
%Mean update
F = eye(size(x, 1)) + (T/2) * Somega(omega);
x = F * x;
%Cov update
G = (T/2) * Sq(x);
P = F * P * F' + G * Rw * G';

end