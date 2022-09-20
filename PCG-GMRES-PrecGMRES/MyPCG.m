function [x, resvec, iter]= mypcg(A,b,tol,maxit,L)

    x(:,1)= zeros(length(b),1);
    r(:,1)= b-A*x(:,1);

    v= L'\r(:,1);
    p(:,1)= L\v;
    v= L'\r(:,1);
    rho(:,1)= r(:,1)'*(L\v);

    k= 1;
    resvec(k)= norm(r(:,k));

    while norm(r(:,k))> tol && k<maxit+1

        z(:,k)= A*p(:,1);
        alpha= rho(:,1)./(((z(:,k)).')*p(:,1));

        x(:,k+1)= x(:,k) + alpha*p(:,1);
        r(:,k+1)= r(:,k)- alpha*z(:,k);

        v= L'\r(:,k+1);
        g(:,k+1)= L\v;
        rho(:,2)= (r(:,k+1).')* g(:,k+1);
        beta= rho(:,2)./rho(:,1);
        p(:,2)= g(:,k+1) + beta*p(:,1);

        p(:,1)= p(:,2);
        rho(:,1)= rho(:,2);

        k= k+1;

        resvec(k)= norm(r(:,k));
    end

    iter= k-1;
end
