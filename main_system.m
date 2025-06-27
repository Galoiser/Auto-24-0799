clear
clc

% System parameters

A = [ 0.9 0 0  0  0.1 0
    0 0.8 0  0.2  0  0
    0 0 0.7 0.1  0  0.1
    0 0.2 0.2  0.7  0  0
    0.1 0 0 0  0.7  0.2
    0 0 0.1 0 0.2 0.7];

C = [1 0 0 0 0 0
     0 1 0 0 0 0
     0 0 1 0 0 0];
 
Q = 0.1*eye(6);
R = 0.11*eye(3);

% Kalman filter parameters
[K,P] = kfilter(A,C,Q,R);


% Problem II and III constraint
delta = 0.25;

% Problem IV weights
w1 = 0.8;
w2 = 1;


save sys.mat


