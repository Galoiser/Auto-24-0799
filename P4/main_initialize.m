% Initialize the endpoints of algorithm 1
clear
clc
load ../sys.mat

% load alarm data
mu_a = [1 2 3]'; % attack mean
VA = diag([0.01 0.1 1]); % attack covariance

% alarm weighting
AK = A-A*K*C;
X = dlyap(AK',A*K*(VA+R)*K'*A'+Q);

 % Problem II:      minimize MAR 
    %             suject to FAR <= delta 
    
% basic parameters
AK = A-A*K*C;
X1 = C*P*C'+R;
X2 = C*X*C'+VA+R;
L1 = (eye(3)+C*inv(AK - eye(6))*A*K)*mu_a;
L = L1*L1';

% SDP
alpha_low = 0;
beta_low = 0;

cvx_begin SDP
variable r(1,1)
minimize(r)
L1*L1'-2*r*X1<=0;
r>=0;
cvx_end

alpha_upp = sqrt(r);

cvx_begin SDP
variable r(1,1)
minimize(r)
L1*L1'-2*r*X2<=0;
r>=0;
cvx_end

beta_upp = sqrt(r);

save endpoints.mat alpha_low alpha_upp beta_low beta_upp 

