function [out]= Alg1(U,V, fun, opts, varargin)

%  min  F(U,V)= -Re(tr(U'AVN))  (U,V) \in St(m,p,C)*St(n,p,C)
%                                U'*U = I_p, where U \in C^{m,p}
%                                V'*V = I_p, where V \in C^{n,p}

%  该程序采用 QR分解作为 Retraction）
%-------------------------------------------------------------------------
% Input:
%           U --- m by p matrix such that U'*U = I
%           V --- n by p matrix such that V'*V=I
%         fun --- objective function and its gradient:
%                 [F, G] = fun(X,  data1, data2)
%                 F, G are the objective function value and gradient, repectively
%                 data1, data2 are addtional data, and can be more
%                 Calling syntax:
%                   [X, out]= OptStiefelGBB(X0, @fun, opts, data1, data2);
%
%        opts --- option structure with fields:
%                 record = 0, no print out
%                 mxitr       max number of iterations
%                 xtol        stop control for ||X_k - X_{k-1}||
%                 gtol        stop control for the projected gradient
%                 ftol        stop control for |F_k - F_{k-1}|/(1+|F_{k-1}|)
%                             usually, max{xtol, gtol} > ftol
% Output:
%           X --- solution
%         Out --- output information
%-------------------------------------------------------------------------
if exist('varargin','var')
    nARG      =    length(varargin);
    if nARG   >    0
        A = varargin{1};
        N = varargin{2};
    end
end



[m,p]=size(U);
[n,p]=size(V);
n2=n/2;
% copy parameters
xtol = opts.xtol;
gtol = opts.gtol;
ftol = opts.ftol;
lambda   = opts.lambda;
delta  = opts.delta;
gamma = 0.85;
alpha = opts.alpha;
record = opts.record;
nt = opts.nt;   crit = ones(nt, 3);


%% Initial function value and gradient
% prepare for iterations
[F,  Gu, Gv] = feval(fun, U, V, varargin{:});  out.nfe = 1;
UG=U'*Gu;    Gradu=Gu-0.5*U*(UG+UG');   Zu=-Gradu;%黎曼梯度
VG=V'*Gv;    Gradv=Gv-0.5*V*(VG+VG');  Zv=-Gradv;%黎曼梯度
nrmGradu=norm(Gradu,'fro');
nrmGradv=norm(Gradv,'fro');
nrmGuv = sqrt(nrmGradu^2+nrmGradv^2);

Q = 1; Cval = F;
%% -------------------------------------------------------------------------
if (opts.record == 1)
    fid = 1;
    fprintf(fid, '%4s %8s %10s %10s %10s %10s %4s\n', 'Iter', 'F(X)', 'nrmGradv', 'nrmGradp',  'xDiff', 'FDiff', 'nls');
end

%% main iteration
for itr = 1 : opts.mxitr
    U_old = U;     V_old = V;  Zu_old=Zu;  Zv_old=Zv;  Gradu_old= Gradu;  Gradv_old= Gradv; F_old=F;
    % scale step size

    nls = 1; deriv = delta*(nrmGradu^2+nrmGradv^2); %deriv
    while 1

        [F,  Gu, Gv] = feval(fun, U, V, varargin{:});
        out.nfe = out.nfe + 1;


        if F <= Cval - alpha*deriv || nls >= 5
            break;
        end
        alpha = lambda*alpha;          nls = nls+1;
    end

    %%  QR分解
    Utan=U_old+alpha*Zu_old;
    Utan1=Utan(1:n2,:);
    Utan2=Utan(n2+1:n,:);
    Utan=Utan1+1i*Utan2;
    [U,Ru]=qr(Utan,0);
    U=[real(U);imag(U)];
    Vtan=V_old+alpha*Zv_old;
    Vtan1=Vtan(1:n2,:);
    Vtan2=Vtan(n2+1:n,:);
    Vtan=Vtan1-1i*Vtan2;
    [V,Rv]=qr(Vtan,0);
    V=[real(V);-imag(V)];
    [F,  Gu, Gv] = feval(fun, U, V, varargin{:});
    UG=U'*Gu;    Gradu=Gu-0.5*U*(UG+UG');   %黎曼梯度
    VG=V'*Gv;    Gradv=Gv-0.5*V*(VG+VG');  %黎曼梯度
    Zu=-Gradu;   Zv=-Gradv;

    Su = U - U_old;   Sv = V - V_old;
    Yu = Zu - Zu_old;    Yv = Zv - Zv_old;
    UVDiff = norm(Su,'fro')/sqrt(m)+norm(Sv,'fro')/sqrt(n);

    FDiff = abs(F_old-F)/(abs(F_old)+1);
    %------------BB步长---加速----------------------------------
    SY = abs(sum(sum(Su.*Yu)))+abs(sum(sum(Sv.*Yv)));
    S=abs(sum(sum(Su.*Su)))+abs(sum(sum(Sv.*Sv)));

    if mod(itr,2)==0
        alpha = sum(sum(S.*S))/SY;
    else
        alpha  = SY/(sum(sum(Yu.*Yu))+sum(sum(Yv.*Yv)));
    end
    alpha = max(min(alpha, 1e20), 1e-20);
    %---------------------------------------------------------
    %
    if (record >= 1)
        fprintf('%4d  %4.3e  %3.2e %3.2e  %3.2e  %3.2e  %2d\n', ...
            itr,  F, nrmGradu, nrmGradv, UVDiff, FDiff, nls);
    end

    crit(itr,:) = [nrmGuv, UVDiff, FDiff];
    mcrit = mean(crit(itr-min(nt,itr)+1:itr, :),1);
    out.feasi_u = norm(U'*U-eye(p),'fro'); out.feasi_v = norm(V'*V-eye(p),'fro');

    if   (UVDiff < xtol && FDiff < ftol)|| nrmGuv < gtol || all(mcrit(2:3) < 10*[xtol, ftol])
        if itr <= 2
            ftol = 0.1*ftol;
            xtol = 0.1*xtol;
            gtol = 0.1*gtol;
        else
            out.msg = 'converge';
            break;
        end
    end
    Qp = Q; Q = gamma*Qp + 1; Cval = (gamma*Qp*Cval + F)/Q;

end

if itr >= opts.mxitr
    out.msg = 'exceed max iteration';
end

out.feasi_u = norm(U'*U-eye(p),'fro'); out.feasi_v = norm(V'*V-eye(p),'fro');
if  out.feasi_u > 1e-13 &&  out.feasi_v > 1e-13
    U = MGramSchmidt(U); V = MGramSchmidt(V);
    [F,~,~] = feval(fun, U, V, varargin{:});
    out.nfe = out.nfe + 1;
    out.feasi_u = norm(U'*U-eye(p),'fro'); out.feasi_v = norm(V'*V-eye(p),'fro');
end

%% output----------
out.itr = itr;
out.U=U;
out.V=V;






