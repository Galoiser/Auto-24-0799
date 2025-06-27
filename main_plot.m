clear
clc
% Compare ROC curves
load sys.mat

% load attack data
mu_a = [1 2 3]'; % attack mean
VA = diag([0.01 0.1 1]); % attack covariance

% alarm weighting
AK = A-A*K*C;
X = dlyap(AK',A*K*(VA+R)*K'*A'+Q);

% Problem I: Data processing for AUC plot
load P1/design1.mat
aw1 = opt_aw;
if aw1'*(eye(3)+C*inv(AK - eye(6))*A*K)*mu_a<0
    aw1 = -aw1;
end
mu_r = 0;
mu_r1 = aw1'*(eye(3)+C*inv(AK - eye(6))*A*K)*mu_a;
sig_r =  sqrt(aw1'*(C*P*C'+R)*aw1);
sig_r1 =  sqrt(aw1'*(C*X*C'+VA+R)*aw1);

FAR_list1 = [];
MAR_list1 = [];
for theta1 = 0:0.01:mu_r1
    FAR1 = 1 - 0.5*(1 + erf((theta1-mu_r)/(sqrt(2)*sig_r)));
    MAR1 = 0.5*(1 + erf((theta1-mu_r1)/(sqrt(2)*sig_r1)));
    FAR_list1 = [FAR_list1 FAR1];
    MAR_list1 = [MAR_list1 MAR1];
end

% Problem II: Data processing for AUC plot
load P2/design2.mat
aw2 = opt_aw;
theta2 = opt_theta;
if aw2'*(eye(3)+C*inv(AK - eye(6))*A*K)*mu_a<0
    aw2 = -aw2;
end
mu_r2 = 0;
mu_r12 = aw2'*(eye(3)+C*inv(AK - eye(6))*A*K)*mu_a;
sig_r2 =  sqrt(aw2'*(C*P*C'+R)*aw2);
sig_r12 =  sqrt(aw2'*(C*X*C'+VA+R)*aw2);

opt_FAR2 = 1 - 0.5*(1 + erf((theta2-mu_r2)/(sqrt(2)*sig_r2)));
opt_MAR2 = 0.5*(1 + erf((theta2-mu_r12)/(sqrt(2)*sig_r12)));

FAR_list2 = [];
MAR_list2 = [];
for theta2 = 0:0.01:mu_r12
    FAR2 = 1 - 0.5*(1 + erf((theta2-mu_r2)/(sqrt(2)*sig_r2)));
    MAR2 = 0.5*(1 + erf((theta2-mu_r12)/(sqrt(2)*sig_r12)));
    FAR_list2 = [FAR_list2 FAR2];
    MAR_list2 = [MAR_list2 MAR2];
end

% Problem III: Data processing for AUC plot
load P3/design3.mat
aw3 = opt_aw;
theta3 = opt_theta;
if aw3'*(eye(3)+C*inv(AK - eye(6))*A*K)*mu_a<0
    aw3 = -aw3;
end
mu_r3 = 0;
mu_r13 = aw3'*(eye(3)+C*inv(AK - eye(6))*A*K)*mu_a;
sig_r3 =  sqrt(aw3'*(C*P*C'+R)*aw3);
sig_r13 =  sqrt(aw3'*(C*X*C'+VA+R)*aw3);

opt_FAR3 = 1 - 0.5*(1 + erf((theta3-mu_r3)/(sqrt(2)*sig_r3)));
opt_MAR3 = 0.5*(1 + erf((theta3-mu_r13)/(sqrt(2)*sig_r13)));

FAR_list3 = [];
MAR_list3 = [];
for theta3 = 0:0.01:mu_r13
    FAR3 = 1 - 0.5*(1 + erf((theta3-mu_r3)/(sqrt(2)*sig_r3)));
    MAR3 = 0.5*(1 + erf((theta3-mu_r13)/(sqrt(2)*sig_r13)));
    FAR_list3 = [FAR_list3 FAR3];
    MAR_list3 = [MAR_list3 MAR3];
end

% Problem IV: Data processing for AUC plot
load P4/design4.mat
aw4 = opt_aw;
theta4 = opt_theta;
if aw4'*(eye(3)+C*inv(AK - eye(6))*A*K)*mu_a<0
    aw4 = -aw4;
end
mu_r4 = 0;
mu_r14 = aw4'*(eye(3)+C*inv(AK - eye(6))*A*K)*mu_a;
sig_r4 =  sqrt(aw4'*(C*P*C'+R)*aw4);
sig_r14 =  sqrt(aw4'*(C*X*C'+VA+R)*aw4);

opt_FAR4 = 1 - 0.5*(1 + erf((theta4-mu_r4)/(sqrt(2)*sig_r4)));
opt_MAR4 = 0.5*(1 + erf((theta4-mu_r14)/(sqrt(2)*sig_r14)));

% tangent line parameters
tang = w1*opt_FAR4+w2*opt_MAR4; % w1*x + w2*y = tang
x_tang = tang/w1;
y_tang = tang/w2;
% 切线 x+y = opt_FAR4+opt_MAR4

FAR_list4 = [];
MAR_list4 = [];
for theta4 = 0:0.01:mu_r14
    FAR4 = 1 - 0.5*(1 + erf((theta4-mu_r4)/(sqrt(2)*sig_r4)));
    MAR4 = 0.5*(1 + erf((theta4-mu_r14)/(sqrt(2)*sig_r14)));
    FAR_list4 = [FAR_list4 FAR4];
    MAR_list4 = [MAR_list4 MAR4];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(FAR_list1,MAR_list1,'LineWidth',1.8)
hold on
plot(FAR_list2,MAR_list2,'LineWidth',1.8)
hold on

plot(FAR_list3,MAR_list3,'LineWidth',1.8)
hold on
plot(FAR_list4,MAR_list4,'LineWidth',1.8)
hold on
scatter(opt_FAR4,opt_MAR4,30,'o','filled','MarkerEdgeColor',"#FF00FF",...
              'MarkerFaceColor',"#FF00FF")
hold on
plot([0, x_tang], [y_tang, 0],'--','Color',[0.15 0.15 0.15]);  % 绘制水平辅助线，'--b' 表示蓝色虚线


scatter(opt_FAR2,opt_MAR2,30,'o','filled')
hold on
scatter(opt_FAR3,opt_MAR3,30,'o','filled')
hold on 

plot([0, opt_FAR2], [opt_MAR2, opt_MAR2], '--','Color',[0.5 0.5 0.5]);  % 绘制水平辅助线，'--b' 表示蓝色虚线
plot([opt_FAR2, opt_FAR2], [0, opt_MAR2], '--','Color',[0.5 0.5 0.5]);  % 绘制垂直辅助线，'--b' 表示蓝色虚线

plot([0, opt_FAR3], [opt_MAR3, opt_MAR3], '--','Color',[0.5 0.5 0.5]);  % 绘制水平辅助线，'--b' 表示蓝色虚线
plot([opt_FAR3, opt_FAR3], [0, opt_MAR3], '--','Color',[0.5 0.5 0.5]);  % 绘制垂直辅助线，'--b' 表示蓝色虚线

legend(sprintf('ROC curve for Problem I (AUC = %.4f)', auc_p1), ...
       sprintf('ROC curve for Problem II (AUC = %.4f)', auc_p2), ...
       sprintf('ROC curve for Problem III (AUC = %.4f)', auc_p3), ...
       sprintf('ROC curve for Problem IV (AUC = %.4f)', auc_p4))
set(gca,'FontSize',14);

xlabel('FAR','Interpreter','latex','fontweight','bold','fontsize',15)
ylabel('MAR','Interpreter','latex','fontweight','bold','fontsize',15)
xlim([0 0.5])
ylim([0 0.5])
hold off

% 图宽6 高5 全屏 解析度300 字体TNR