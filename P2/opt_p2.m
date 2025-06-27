% Optimization procedure corresponding to Theorem 3 (bisection method)

function [w_dagger,theta_dagger,check] = opt_p2(A,C,Q,R,X,K,P,mu_a,VA,delta)
    AK = A-A*K*C;
    X1 = C*P*C'+R;
    X2 = C*X*C'+VA+R;
    L1 = (eye(3)+C*inv(AK - eye(6))*A*K)*mu_a;
    L = L1*L1';
    
    eta1 = 0;
    eta2 = 20;

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
        sqrt(2)*eta3*a+sqrt(2)*erfinv(1-2*delta)*b<=r;
        cvx_end
        if r > 1
            eta2 = eta3;
        else
            eta1 = eta3;
        end
    end

    if norm(r-1)>1e-3
        error('The FAR constraint <= delta is unachievable!')
    end
    check = r;
    w_dagger = w;
    theta_dagger = sqrt(2)*erfinv(1-2*delta)*norm(sqrtm(X1)*w);
end
