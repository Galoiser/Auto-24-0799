% Main procedure for solving Problem I

clear
clc
load ../sys.mat

% load attack data
mu_a = [1 2 3]'; % attack mean
VA = diag([0.01 0.1 1]); % attack covariance

% alarm weighting
AK = A-A*K*C;
X = dlyap(AK',A*K*(VA+R)*K'*A'+Q);

aw = 0.5-rand(3,1);
aw = aw/norm(aw);
[opt_aw,check] = opt_p1(A,C,Q,R,X,K,P,mu_a,VA);

% AUC Calculation
mu_r = 0;
mu_r1 = aw'*(eye(3)+C*inv(AK - eye(6))*A*K)*mu_a;
sig_r = aw'*(C*P*C'+R)*aw;
sig_r1 = aw'*(C*X*C'+VA+R)*aw;


opt_mu_r = 0;
opt_mu_r1 = opt_aw'*(eye(3)+C*inv(AK - eye(6))*A*K)*mu_a;
opt_sig_r = opt_aw'*(C*P*C'+R)*opt_aw;
opt_sig_r1 = opt_aw'*(C*X*C'+VA+R)*opt_aw;

if opt_mu_r1 < 0
    opt_aw = -opt_aw;
    opt_mu_r = 0;
    opt_mu_r1 = opt_aw'*(eye(3)+C*inv(AK - eye(6))*A*K)*mu_a;
    opt_sig_r = opt_aw'*(C*P*C'+R)*opt_aw;
    opt_sig_r1 = opt_aw'*(C*X*C'+VA+R)*opt_aw;
end

% Results display
% Problem I is threshold-independent, so FAR and MAR are not calculated here.

% Results display
disp('Alarm weight and threshold')
alarm_weight = opt_aw

disp('AUC value')
auc_p1 = 1-normcdf(abs(opt_mu_r-opt_mu_r1)/sqrt(opt_sig_r+opt_sig_r1))

save design1.mat opt_aw auc_p1
