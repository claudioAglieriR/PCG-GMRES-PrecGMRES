function[x, iter, resvec, flag]= myprecgmres(A, b, tol, maxit, x0, L, U)

    [n,n]= size(A);
    flag= 0;
    k= 1;

    x= x0;
    M= L*U;
    r0= M\(b-(A*x(:,k)));

    rho= norm(r0, 2);
    initial_rho= rho;
    beta= norm(rho, 2);
    v(:,k)= r0/beta;
    stopwhile= tol*(norm((L*U)\b,2));

    while (rho>stopwhile) && (k<maxit)

        z= A*v(:,k);
        v(:,k+1)= M\z;

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
        resvec(k)= rho;
        rho= norm(resvec(k),2);
        k= k+1;
    end

    x= x0 + v(:,1:end-1)*z;
    iter= k-1;

    % Computation of the "true" final residual norm

    true_residual= norm((b-A*x),2);
    resvec= [initial_rho, resvec];

    Sentence1= ['\n\n The final computed residual is %9.3e . The "true" final residual norm is %9.3e .'];

    fprintf(Sentence1, (resvec(end)), true_residual)

    if (resvec(end)==0)

        flag= -1;
    end
end
