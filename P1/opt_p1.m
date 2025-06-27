% Optimization procedure corresponding to Theorem 2

function [w_dagger,check] = opt_p1(A,C,Q,R,X,K,P,mu_a,VA)
    AK = A-A*K*C;
    X1 = C*P*C'+R;
    X2 = C*X*C'+VA+R;
    L1 = (eye(3)+C*inv(AK - eye(6))*A*K)*mu_a;
    L = L1*L1';

    cvx_begin SDP
    variable r(1,1) 
    minimize(r)
    L1*L1'-r*(X1+X2)<=0;
    r>=0;
    cvx_end

    [V,D] = eig(L1*L1'-r*(X1+X2));
    w_dagger = 0;
    for i = 1:3
        if norm(D(i,i))<1e-4
            w_dagger = V(:,i);
        end
    end
    check = r;
    w_dagger = w_dagger/norm(w_dagger);

end
