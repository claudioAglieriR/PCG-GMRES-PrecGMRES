clear all


%% Test of my implementation of GMRES (mygmres)


A= readmatrix('mat13041.rig.txt');
A= spconvert(A);
[m,n]= size(A);

x0= zeros(n,1);
xg= zeros(n,1);

itmax= 550;
tol= 1e-10;

for i= 1:1:n

    xg(i,1)= 1/(sqrt(i));
end

b= A*xg;

[myx, myiter, myresvec, myflag]= mygmres(A,b,tol,itmax,x0);

Sentence= ['\n Using GMRES with my implementation, the value of iter is %4.2f and the last value of resvec is %9.3e .' ...
    '   \n The value of flag is %.15g .' ...
    '\n\n\n A plot of the residual norm vs the number of iterations will be shown on Figure 1.'];
fprintf(Sentence, myiter, (myresvec(end)), myflag)

% Plot of the residual norm vs the number of iterations

axes('YScale', 'log')
xlabel('Iteration');
ylabel('Residual');

figure(1)
hold on
semilogy(0:myiter, myresvec,'k--*');
legend('mygmres');
xlabel('Iteration');
ylabel('Residual norm');
hold off
