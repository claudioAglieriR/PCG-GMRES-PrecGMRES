clear all


%% Test of my implementation of Preconditioned GMRES (myprecgmres)


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

% Computation of the ILU preconditioner

setup.type= 'crout';
setup.droptol= 0.1;
[L,U]= ilu(A, setup);

[myx, myiter, myresvec, myflag]= myprecgmres(A,b,tol,itmax,x0,L,U);

% Exercise 6.(c)

[x,flag,relres,iter,resvec]= gmres(A,b,10000,tol,itmax,L,U);

tot_iter=iter(2);

Sentence2= ['\n\n On Figure 1 it is possible to see the plot of my implementation of preconditioned gmres' ...
    ' and of MATLAB gmres.' ...
    '\n\n Figure 2 is a table showing the results obtained with my implementation of preconditioned gmres' ...
    ' and with MATLAB gmres.'];

fprintf(Sentence2)

axes('YScale', 'log')
xlabel('Iteration');
ylabel('Residual norm');

figure(1)

hold on
semilogy(0:tot_iter, resvec,'r-*', 0:myiter, myresvec,'g-o');
legend('gmres','myprecgmres');
xlabel('Iteration ');
ylabel('Residual norm');
hold off

% Exercise 6.(d)

figure(2)

hold on
Method= {'myprecgmres'; 'gmres'};
iter= [myiter; tot_iter];
resvec= [myresvec(end); resvec(end)];
flag= [myflag; flag];
T = table(iter, resvec, flag, 'RowNames',Method);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

hold off
