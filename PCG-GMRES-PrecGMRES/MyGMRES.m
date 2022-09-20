function[x, iter, resvec, flag]= mygmres(A, b, tol, maxit, x0)

    [n,n]= size(A);
    flag= 0;
    k= 1;
    
    x(:,k)= x0;
    r0= b-(A*x(:,k));
    rho= norm(r0, 2);
    resvec(k)= rho;
    beta= rho;
    v(:,k)= r0./beta;

    stopwhile= (tol*norm(b,2));

    while (rho>stopwhile) && (k<maxit)

        v(:,k+1)= A*v(:,k);

        for j=1:1:k

            H(j, k)= (v(:,k+1).')*v(:,j);
            v(:,k+1)= v(:,k+1)- (H(j, k))*v(:,j);
        end

        H(k+1, k)= norm(v(:,k+1),2);
        v(:, k+1)= v(:, k+1)./H(k+1, k);


        [m1, n1]= size(H);
        MAT= speye(m1);
        e_1= MAT(:,1);
        z= lsqminnorm(H,beta.*e_1);
        rho= norm(beta.*e_1-H*z);

        k= k+1;

        % resvec is the vector with norm of the residuals

        resvec(k)= rho;
    end

    % x is the solution vector

    x= x0 + v(:,1:end-1)*z;

    % iter is the number of iterations employed

    iter= k-1;

    if (resvec(end)==0)

        flag= -1;
    end
end
