% Find the optimal alpha and beta (Algorithm 1)

function [w_dagger,theta_dagger,alpha,beta,obj,J_grad] = opt_p4(A,C,Q,R,X,K,P,mu_a,VA,w1,w2,alpha_low,alpha_upp,beta_low,beta_upp)
    AK = A-A*K*C;
    X1 = C*P*C'+R;
    X2 = C*X*C'+VA+R;
    L1 = (eye(3)+C*inv(AK - eye(6))*A*K)*mu_a;
    L = L1*L1';
    
    eta1 = alpha_low;
    eta2 = alpha_upp;
    J_grad = 1;


    while norm(J_grad)>1e-4
        eta3 = (eta1+eta2)/2;
        alpha = eta3;
        [w_dagger,theta_dagger,check,beta,obj] = opt_beta(A,C,Q,R,X,K,P,mu_a,VA,alpha,w1,w2,beta_low,beta_upp);
        h_grad = norm(sqrtm(X1)*w_dagger)/norm(sqrtm(X2)*w_dagger);
        J_grad = (w1/2)*erf_grad(alpha)-(w2/2)*erf_grad(beta)*h_grad;
        if J_grad >0
            eta1 = eta3;
        else
            eta2 = eta3;
        end
    end
end
