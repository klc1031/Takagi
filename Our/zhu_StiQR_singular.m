function [out]= zhu_StiQR_singular(U, fun, opts, varargin)
% Riemannian nonmonotone conjugate gradient method for Problem
%-------------------------------------------------------------------------
% problem model：
%  min -1/2 tr(Jp U'*Jn*A*U*N) ,
%  s.t.   U'*U = I2p

% where U is a 2n-by-2p matrix, Jn = [0 In; In 0],
% I2p is the 2p-by-2p identity matrix.


%-------------------------------------------------------------------------
% Input:
%           U --- n by p matrix such that U'*U = I2p
%         fun --- objective function and its gradient:
%                 [F, G] = fun(X,  data1, data2)
%                 F, G are the objective function value and gradient, repectively
%                 data1, data2 are addtional data, and can be more
%                 Calling syntax:
%                   [out]= zhu_StiQR_singular(U0, @fun, opts, data1, data2);
%
%        opts --- option structure with fields:
%                 record = 0, no print out
%                 mxitr       max number of iterations
%                 xtol        stop control for ||U_k - U_{k-1}||
%                 gtol        stop control for the projected gradient
%                 ftol        stop control for |F_k - F_{k-1}|/(1+|F_{k-1}|)
%                             usually, max{xtol, gtol} > ftol
%    varargin --- parameter
%                 A         2n-by-2n matrix
%                 N         truncated matrix
%                 J1        permutation matrix Jn
%                 J2        permutation matrix Jp
%

% Output:
%         Out --- output information
%-------------------------------------------------------------------------


[n,p]=size(U);
n2=n/2;
p2=p/2;
In=eye(n2);
% copy parameters
xtol = opts.xtol;
gtol = opts.gtol;
ftol = opts.ftol;
delta  = opts.delta;
lambda   = opts.lambda;
mm=opts.mm;
iscomplex=1;
alpha = opts.alpha;
record = opts.record;
nt = opts.nt;   crit = ones(nt, 3);

%% Initial function value and gradient
% prepare for iterations
[F, Gu] = feval(fun, U, varargin{:});  out.nfe = 1;
UG=U'*Gu;  Gradu=Gu-0.5*U*(UG+UG');% Riemannian gradient

Zu=-Gradu; %search direction
prodGZ=sum(dot(Gradu,Zu,1));
nrmGradu=norm(Gradu,'fro');

Fm=F;


%-------------------------------------------------------------------------
if (opts.record == 1)
    fid = 1;
    fprintf(fid, '%4s %8s %10s %10s %10s %4s\n', 'Iter', 'F(X)', 'nrmGu',  'XDiff', 'FDiff', 'nls');
end

%% main iteration
for itr = 1 : opts.mxitr
    U_old=U;   F_old = F;
    nrmGradu_old = nrmGradu;
    Gradu_old= Gradu;
    Zu_old=Zu;

    % scale step size
    nls = 1;
    while 1

        %  QR-based retraction
        Utan=U_old+alpha*Zu_old;
        Utan1=Utan(1:n2,1:p2);
        Utan2=Utan(1:n2,p2+1:p);
        Utan=Utan1+1i*Utan2;
        [U1,Ru]=qr(Utan,0);
        U=[real(U1),imag(U1);-imag(U1),real(U1)];


        [F,Gu] = feval(fun, U, varargin{:});
        out.nfe = out.nfe + 1;

        if F <= max(Fm) + alpha*delta*prodGZ || nls >= 5
            break;
        end
        alpha = lambda*alpha;          nls = nls+1;

    end

    % vector transport
    Zuu=norm(Zu_old,'fro');
    Zu_old1=Zu_old(1:n2,1:p2);
    Zu_old2=Zu_old(1:n2,p2+1:p);
    Zu_old3=Zu_old1+1i*Zu_old2;
    X=U1'*Zu_old3*inv(Ru);
    skewh=tril(X,-1)-tril(X,-1)'+1i*diag(diag(X));
    Tzu1=U1*skewh+(In-U1*U1')*Zu_old3*inv(Ru);
    TZu=[real(Tzu1),imag(Tzu1);-imag(Tzu1),real(Tzu1)];
    TZuu=norm(TZu,'fro');

    if TZuu >  Zuu
        TZu=Zuu/TZuu*TZu; % scaling
    end


    UG=U'*Gu;  Gradu=Gu-0.5*U*(UG+UG');
    nrmGradu=norm(Gradu,'fro');
    nrmGu = sqrt(nrmGradu^2);

    mm1=sum(dot(Gradu,TZu,1));

    betaD = (nrmGradu^2)/max([mm1-prodGZ, -prodGZ]);
    betaFR = (nrmGradu^2)/(nrmGradu_old^2);
    beta = min(betaD,betaFR);%beta

    Zu=-Gradu+beta*TZu;
    prodGZ=sum(dot(Gradu,Zu,1));

    Su = U - U_old;
    UDiff = norm(Su,'fro')/sqrt(n);
    FDiff = abs(F-F_old)/(abs(F_old)+1);

    if (record >= 1)
        fprintf('%4d  %4.3e  %3.2e  %3.2e  %3.2e  %2d\n', ...
            itr,  F, nrmGradu, UDiff, FDiff, nls);
    end

    crit(itr,:) = [nrmGu, UDiff, FDiff];
    mcrit = mean(crit(itr-min(nt,itr)+1:itr, :),1);

    if ( UDiff < xtol && FDiff < ftol ) || nrmGu < gtol || all(mcrit(2:3) < 10*[xtol, ftol])
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

    TT(itr)=nrmGu;
    FF(itr)=F;

end

if itr >= opts.mxitr
    out.msg = 'exceed max iteration';
end

out.feasi_u = norm(U'*U-eye(p),'fro');
if  out.feasi_u > 1e-13
    U = MGramSchmidt(U);
    [F,~] = feval(fun, U, varargin{:});
    out.nfe = out.nfe + 1;
    out.feasi_u = norm(U'*U-eye(p),'fro');
end

%% output----------
out.gradfu = nrmGradu;
out.fval = F;
out.Gradu=nrmGu;
out.itr = itr;
out.TT=TT;
out.FF=FF;
out.U=U;



