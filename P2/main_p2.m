% Main procedure for solving Problem II

clear
clc

% main problem2 
load ../sys.mat

% load attack data
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


[w_dagger,theta_dagger,check] = opt_p2(A,C,Q,R,X,K,P,mu_a,VA,delta);

% scaling parameter
lambda = 1/norm(w_dagger);

opt_aw = lambda*w_dagger;
opt_theta = lambda*theta_dagger;

% FAR and MAR calculation
opt_mu_r = 0;
opt_mu_r1 = opt_aw'*(eye(3)+C*inv(AK - eye(6))*A*K)*mu_a;
opt_sig_r = sqrt(opt_aw'*(C*P*C'+R)*opt_aw);
opt_sig_r1 = sqrt(opt_aw'*(C*X*C'+VA+R)*opt_aw);

% Results display
disp('Alarm weight and threshold')
alarm_weight = opt_aw
alarm_threshold = opt_theta

disp('FAR and MAR values')
FAR = 1 - 0.5*(1 + erf((opt_theta-opt_mu_r)/(sqrt(2)*opt_sig_r)))
MAR = 0.5*(1 + erf((opt_theta-opt_mu_r1)/(sqrt(2)*opt_sig_r1)))

disp('weighte sum of FAR and MAR')
weighted_sum = w1*FAR+w2*MAR

disp('AUC value')
auc_p2 = 1-normcdf(abs(opt_mu_r-opt_mu_r1)/sqrt(opt_sig_r^2+opt_sig_r1^2))

save design2.mat opt_aw opt_theta auc_p2
