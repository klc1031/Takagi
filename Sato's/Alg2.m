function [out]= Alg2(U,V, fun, opts, varargin)
%-------------------------------------------------------------------------
% Input:
%           V --- n by p matrix such that V'*V = I
%           P --- 2P by p matrix such that P'*P=I
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
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------
[m,p]=size(U);
[n,p]=size(V);
% copy parameters
xtol = opts.xtol;
gtol = opts.gtol;
ftol = opts.ftol;
delta  = opts.delta;
lambda   = opts.lambda;
mm=opts.mm;
alpha = opts.alpha;
record = opts.record;
nt = opts.nt;   crit = ones(nt, 3);

%% Initial function value and gradient
% prepare for iterations
[F,  Gu, Gv] = feval(fun, U, V, varargin{:});  out.nfe = 1;
UG=U'*Gu;  Gradu=Gu-0.5*U*(UG+UG');%%%%ÀèÂüÌÝ¶È
VG=V'*Gv;  Gradv=Gv-0.5*V*(VG+VG');


Zu=-Gradu; Zv=-Gradv;%%ËÑË÷·½Ïò
prodGZ=sum(dot(Gradu,Zu,1))+sum(dot(Gradv,Zv,1));
nrmGradu=norm(Gradu,'fro'); nrmGradv=norm(Gradv,'fro');

Fm=F;


%-------------------------------------------------------------------------
if (opts.record == 1)
    fid = 1;
    fprintf(fid, '%4s %8s %10s %10s %10s %10s %4s\n', 'Iter', 'F(X)', 'nrmGu', 'nrmGv',  'XDiff', 'FDiff', 'nls');
end

%% main iteration
for itr = 1 : opts.mxitr
    U_old = U;  V_old=V;   F_old = F;
    nrmGradu_old = nrmGradu;  nrmGradv_old = nrmGradv;
    Gradu_old= Gradu;  Gradv_old= Gradv;
    Zu_old=Zu;  Zv_old=Zv;
    % scale step size

    nls = 1;
    while 1

        Utan=U_old+alpha*Zu_old;  [U,Ru]=qr(Utan,0);
        Vtan=V_old+alpha*Zv_old;  [V,Rv]=qr(Vtan,0);

        [F,Gu,Gv] = feval(fun, U, V, varargin{:});
        out.nfe = out.nfe + 1;

        if F <= max(Fm) + alpha*delta*prodGZ || nls >= 5
            break;
        end
        alpha = lambda*alpha;          nls = nls+1;

    end
    ZRinu=Zu/Ru;
    UZRinu=U'*ZRinu;
    UUZRinu=U*UZRinu;
    Uu=tril(UZRinu,-1);
    TZu=U*(Uu-Uu')+ZRinu-UUZRinu;

    ZRinv=Zv/Rv;
    VZRinv=V'*ZRinv;
    VVZRinv=V*VZRinv;
    Uv=tril(VZRinv,-1);
    TZv=V*(Uv-Uv')+ZRinv-VVZRinv;


    UG=U'*Gu;  Gradu=Gu-0.5*U*(UG+UG');
    VG=V'*Gv;  Gradv=Gv-0.5*V*(VG+VG');
    nrmGradu=norm(Gradu,'fro'); nrmGradv=norm(Gradv,'fro');
    nrmGuv = sqrt(nrmGradu^2+nrmGradv^2);

    mm1=sum(dot(Gradu,TZu,1))+sum(dot(Gradv,TZv,1));

    betaD = (nrmGradu^2+nrmGradv^2)/max([mm1-prodGZ, -prodGZ]);
    betaFR = (nrmGradu^2+nrmGradv^2)/(nrmGradu_old^2+nrmGradv_old^2);
    beta = min(betaD,betaFR);
    %

    Zu=-Gradu+beta*TZu; Zv=-Gradv+beta*TZv;
    prodGZ=sum(dot(Gradu,Zu,1))+sum(dot(Gradv,Zv,1));

    Su = U - U_old;   Sv = V - V_old;
    UVDiff = norm(Su,'fro')/sqrt(m)+norm(Sv,'fro')/sqrt(n);
    FDiff = abs(F-F_old)/(abs(F_old)+1);

    Su=alpha*(Zu_old); Sv=alpha*(Zv_old);
    Yu=Gradu-Gradu_old; Yv=Gradv-Gradv_old;




    if (record >= 1)
        fprintf('%4d  %4.3e  %3.2e %3.2e  %3.2e  %3.2e  %2d\n', ...
            itr,  F, nrmGradu, nrmGradv, UVDiff, FDiff, nls);
    end

    crit(itr,:) = [nrmGuv, UVDiff, FDiff];
    mcrit = mean(crit(itr-min(nt,itr)+1:itr, :),1);
    if ( UVDiff < xtol && FDiff < ftol ) || nrmGuv < gtol || all(mcrit(2:3) < 10*[xtol, ftol])
        if itr <= 2
            ftol = 0.1*ftol;
            xtol = 0.1*xtol;
            gtol = 0.1*gtol;
        else
            out.msg = 'converge';
            break;
        end
    end

    Q=Fm;
    Fm=[F,Q];
    mk=min(mm-1,itr);
    Fm=Fm(1:mk+1);

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



