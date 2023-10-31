function [out]= wen_Stiefel_singular(U, fun, opts, varargin)
% Steepest descent method for Problem
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
%                   [out]= wen_Stiefel_singular(U0, @fun, opts, data1, data2);
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
% Output:
%         Out --- output information
%-------------------------------------------------------------------------



if exist('varargin','var')
    nARG      =    length(varargin);
    if nARG   >    0
        A = varargin{1};
        N = varargin{2};
        J1 = varargin{3};
        J2 = varargin{4};
    end
end



[n,p]=size(U);
n2=n/2;
p2=p/2;
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
[F, Gu] = feval(fun, U, varargin{:});  out.nfe = 1;
UG=U'*Gu;    Gradu=Gu-0.5*U*(UG+UG');  %Riemannian gradient
Zu=-Gradu;
nrmGradu=norm(Gradu,'fro');


Q = 1; Cval = F;
%% -------------------------------------------------------------------------
if (opts.record == 1)
    fid = 1;
    fprintf(fid, '%4s %8s %10s %10s %10s %10s %4s\n', 'Iter', 'F(X)', 'nrmGradv', 'nrmGradp',  'xDiff', 'FDiff', 'nls');
end

%% main iteration
for itr = 1 : opts.mxitr
    U_old = U;   Zu_old=Zu;   Gradu_old= Gradu; F_old=F;

    % scale step size
    nls = 1; deriu = delta*(nrmGradu^2);
    while 1

        [F,  Gu] = feval(fun, U, varargin{:});
        out.nfe = out.nfe + 1;


        if F <= Cval - alpha*deriu || nls >= 5
            break;
        end
        alpha = lambda*alpha;          nls = nls+1;
    end

    %  QR-based retraction
    Utan=U_old+alpha*Zu_old;
    Utan1=Utan(1:n2,1:p2);
    Utan2=Utan(1:n2,p2+1:p);
    Utan=Utan1+1i*Utan2;
    [U,Ru]=qr(Utan,0);

    U=[real(U),imag(U);-imag(U),real(U)];
    [F,  Gu] = feval(fun,  U, varargin{:});
    UG=U'*Gu;    Gradu=Gu-0.5*U*(UG+UG');
    Zu=-Gradu;

    nrmGradu=norm(Gradu,'fro');


    Su = U - U_old;
    Yu = Zu - Zu_old;
    UDiff = norm(Su,'fro')/sqrt(n);

    FDiff = abs(F_old-F)/(abs(F_old)+1);

    if (record >= 1)
        fprintf('%4d  %4.3e   %3.2e  %3.2e  %3.2e  %2d\n', ...
            itr,  F,  nrmGradu, UDiff, FDiff, nls);
    end

    crit(itr,:) = [nrmGradu, UDiff, FDiff];
    mcrit = mean(crit(itr-min(nt,itr)+1:itr, :),1);
    out.feasi_u = norm(U'*U-eye(p),'fro');

    if   UDiff < xtol || FDiff < ftol|| all(mcrit(2:3) < 10*[xtol, ftol])
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


    TT(itr)=nrmGradu;
    FF(itr)=F;

end

if itr >= opts.mxitr
    out.msg = 'exceed max iteration';
end

out.feasi_u = norm(U'*U-eye(p),'fro');
if   out.feasi_u > 1e-13
    U = MGramSchmidt(U);
    [F,~] = feval(fun,  U, varargin{:});
    out.nfe = out.nfe + 1;
    out.feasi_u = norm(U'*U-eye(p),'fro');
end

%% output----------
out.gradfu = nrmGradu;
out.fval = F;
out.itr = itr;
out.TT=TT;
out.FF=FF;
out.U=U;






