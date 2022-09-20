clear all


%% Test of my implementation of PCG (mypcg)


A= delsq(numgrid('S',102));
L= ichol(A);
n= size(A,1);
b= A*ones(n,1);
tol= 1.e-8;
maxit= 50;

[x, flag, relres, iter, resvec]= pcg(A, b, tol, maxit, L, L');

Sentence1= [' Using PCG with Matlab implementation, the value of iter is %4.2f. The last value of resvec is %9.3e .'];

fprintf(Sentence1, iter, (resvec(end)))


[my_x, my_resvec, my_iter]= mypcg(A, b, tol, maxit, L);

Sentence2= ['\n\n Using PCG with my implementation, the value of iter is %4.2f. The last value of resvec is %9.3e .' ...
    '\n\n\n A plot which compares my PCG with Matlab implementation will be shown on Figure 1.'];
fprintf(Sentence2, my_iter, (my_resvec(end)))


semilogy(0:iter, resvec, 'k--o', 0:my_iter, my_resvec, 'g--*');
legend('PCG','mypcg');
xlabel('Iteration');
ylabel('Residual norm');
