% Inner bisection (Algorithm 1)

function [w_dagger,theta_dagger,check,beta,obj] = opt_beta(A,C,Q,R,X,K,P,mu_a,VA,alpha,w1,w2,beta_low,beta_upp)
    AK = A-A*K*C;
    X1 = C*P*C'+R;
    X2 = C*X*C'+VA+R;
    L1 = (eye(3)+C*inv(AK - eye(6))*A*K)*mu_a;
    L = L1*L1';
    
    eta1 = beta_low;
    eta2 = beta_upp;

    while norm(eta1 - eta2)>1e-6
        eta3 = (eta1+eta2)/2;
        cvx_begin 
        variable r(1,1) 
        variable a(1,1) 
        variable b(1,1) 
        variable w(3,1) 
        minimize(r)
        L1'*w==1;
        norm(sqrtm(X2)*w)<=a;
        norm(sqrtm(X1)*w)<=b;
        sqrt(2)*alpha*b+sqrt(2)*eta3*a<=r;
        cvx_end
        if r > 1
            eta2 = eta3;
        else
            eta1 = eta3;
        end
    end

    if norm(r-1)>1e-3
        error('Problem with initialization!!!!!!')
    end
    check = r;
    w_dagger = w/norm(w);
    theta_dagger =  sqrt(2)*alpha*norm(sqrtm(X1)*w)/norm(w);
    beta = eta3;
    
    obj = (w1/2)*erf(alpha)+(w2/2)*erf(beta);
end
