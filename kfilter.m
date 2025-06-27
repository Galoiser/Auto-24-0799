% Return parameters of the Kalman filter

function [L,P] = kfilter(A,C,Q,R)
P_ = zeros(size(A));
P = ones(size(A));
while norm(P - P_)>1e-5
    P = P_;
    Ph = A*P*A'+Q;
    L = Ph*C'*inv(C*Ph*C'+R);
    P_ = (eye(size(A)) - L*C)*Ph;
end

