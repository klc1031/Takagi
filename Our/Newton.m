function [U,G,iter,Fu,Grad] = Newton(U0, fun,opts,varargin)

if exist('varargin','var')
    nARG      =    length(varargin);
    if nARG   >    0
        A = varargin{1};
        N = varargin{2};
        Jn = varargin{3};
        Jp = varargin{4};
    end
end

%设置初始参数
k_max=50; %最大迭代次数
[n,p]=size(U0);
n2=n/2;
p2=p/2;
xtol = opts.xtol;
gtol = opts.gtol;
ftol = opts.ftol;
U=U0;
[Fu, Gu] = feval(fun, U, varargin{:});
UG=U'*Gu;
Gradu=Gu-0.5*U*(UG+UG');
Grad=norm(Gradu);
F(1)=Fu;
G(1)=Grad;
for k = 1:k_max
    U_old = U;F_old=Fu;
    g=reshape(-Gradu,[n*p,1]);%vec
    [M]=xishu(n,p,U_old,A,N,Jn,Jp);%solve the Newton equation
    xi=M\g;
    xi=reshape(xi,[n,p]);

    Utan=U_old+xi;
    Utan1=Utan(1:n2,1:p2);
    Utan2=Utan(1:n2,p2+1:p);
    Utan=Utan1+1i*Utan2;
    [U1,~]=qr(Utan,0);%retraction
    U=[real(U1),imag(U1);-imag(U1),real(U1)];


    [Fu, Gu] = feval(fun, U, varargin{:});
    UG=U'*Gu;
    Gradu=Gu-0.5*U*(UG+UG');
    Grad=norm(Gradu,'fro');
    Su = U - U_old;
    UDiff = norm(Su,'fro')/sqrt(n);
    FDiff = abs(F_old-Fu)/(abs(F_old)+1);
    if   (UDiff < xtol && FDiff < ftol)|| Grad < gtol
        if k <= 2
            ftol = 0.1*ftol;
            xtol = 0.1*xtol;
            gtol = 0.1*gtol;
        else
            out.msg = 'converge';
            break;
        end
    end
    F(k+1)=Fu;
    G(k+1)=Grad;
end

iter=k;
end