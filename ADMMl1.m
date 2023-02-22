function [M, cnt, fail] = ADMMl1(A, r, mu, Nmax, delta)

[m, n] = size(A);

E_old = zeros(m,n);

Lambda = zeros(m,n);

fail = false;

cnt = 0;

while 1

    cnt = cnt + 1;

    Y_old = E_old - Lambda/mu + A;

    [U_t,S_t,V_t] = svds(Y_old,r);
    U_new = U_t*(S_t^(1/2));
    V_new = (S_t^(1/2))*V_t';

    Y_old = U_new*V_new + Lambda/mu - A;

    E_new = zeros(m,n);

    for i = 1:m
        for j = 1:n
            E_new(i,j) = (max(abs(Y_old(i,j)) - 1/mu,0)/(max(abs(Y_old(i,j)) - 1/mu,0)+1/mu))*Y_old(i,j);
        end
    end

    Lambda = Lambda + mu*(U_new*V_new - E_new - A);

    if (cnt >= Nmax) || ((norm(U_new*V_new - E_new - A, 'fro')/norm(A,'fro'))<delta)
        if (cnt >= Nmax)
            fail = true;
        end
        break
    end

    E_old = E_new;

end

M = U_new*V_new;


end
