function [x,cell] = CR(A,c,itr)
%CR   Conjugate residual method
%   X = CR(A,B)  attempts to solve the system of linear equations A*x=b
%   for X. The n-by-n coefficient matrix A must be symmetric
%    and the right hand side column vector b must have length n.

% Check matrix and right hand side vector inputs have appropriate sizes
[n,~] = size(A);
maxit=5000;
tol_x=1e-6;
tol_f=1e-6;
x = zeros(n,1);
% Initial guess
r=c-A*x;
p=r;
Plot_CR=[];
for k = 1:maxit
    r_old=r;  Ar_old=A*r_old;  Ap=A*p; x_old=x;

    alafa=r_old'*Ar_old/(Ap'*Ap);

    x=x+alafa*p;

    r=r_old-alafa*Ap;

    Ar=A'*r;

    beta=r'*Ar/(r_old'*Ar_old);

    p=r+beta*p;

    xd = x - x_old;
    xDiff = norm(xd,'fro')/norm(x_old,'fro');
    normr=norm(r)/norm(c);

    Plot_CR(k)=normr;

    % Break if update increment is below tolerance
    if  normr< tol_f || xDiff< tol_x
        break;
    end
end
cell{itr}=Plot_CR;