function  [U,V,out]=Alg3(U,V,fun, opts, varargin)
%-------------------------------------------------------------------------

record = opts.record;
[m,p]=size(U);
[n,p]=size(V);
pp=p^2;
eye1 = eye(p);
eye2 = eye(pp);
eye3 = eye(m-p);
eye4 = eye(n-p);

if p==1
    K=zeros(1,1);
else
    K=compute_Kn(p);
end

S=compute_Sn(p);


eye2p = eye(2*p);

%---Used when Cayley transform as Retraction----------
invU = true;
if p < n/2;
    invU = false;
end
%-----------------------------------------------------
if exist('varargin','var')
    nARG      =    length(varargin);
    if nARG   >    0
        A = varargin{1};
        N = varargin{2};
    end
end
%AA=A'*A; AB=A'*B; BA=B'*A; BB=B'*B;
%% -------------------------------------------------------------------------
if (opts.record == 1)
    fid = 1;
    fprintf(fid, '%4s %8s %10s %10s \n', 'Iter', 'F(X)', 'nrmGradv', 'nrmGradp');
end


%% main iteration
for itr = 1 : 100

    [U0, ~] = qr(U);
    Ubot = U0(:, p+1:end);  % generate U_{bot}
    [V0, ~] = qr(V);
    Vbot = V0(:, p+1:end);  % generate V_{bot}

    S1=0.5*(U'*A*V*N+N'*V'*A'*U);   % generate S1
    S2=0.5*(V'*A'*U*N+N'*U'*A*V);   % generate S1



    Tpp=vecperm(p,p);
    %T2=vecperm(p,n-p);

    KK1=0.5*K'*(eye2-Tpp);
    KK2=0.5*S'*(eye2+Tpp);

    W11=kron(S1,eye1);
    WK=W11*K;    WS=W11*S;
    W=[WK,1i*WS];
    %     RW=real(W);
    %     IW=imag(W);

    UAV=U'*A*V;
    Q11=kron(N,UAV);
    QK=Q11*K;  QS=Q11*S;
    Q=[QK,1i*QS];

    UAVt=U'*A*Vbot;
    T11=kron(N,UAVt);
    T=[T11,1i*T11];

    qfx=0.5*p*(p-1); qfy=2*(m-p)*p;  qfx2=0.5*p*(p+1);
    QE1=KK1*real(W);    QE2=KK2*imag(W);
    if p==1
        QF1=zeros(p,2*(m-p));     QF2=zeros(p,2*(m-p));
    else
        QF1=zeros(qfx,qfy);     QF2=zeros(qfx2,qfy);
    end
    QM1=-KK1*real(Q);   QM2=-KK2*imag(Q);
    QN1=-KK1*real(T);   QN2=-KK2*imag(T);


    P11=kron(S1,eye3);
    P=[P11,1i*P11];

    UtAV=Ubot'*A*V;
    KK11=kron(N,UtAV);
    KKK=KK11*K;  KKS=KK11*S;
    KK=[KKK,1i*KKS];

    UtAVt=Ubot'*A*Vbot;
    B11=kron(N,UtAVt);
    B=[B11,1i*B11];

    rex=p*(m-p);
    if p==1
        RE1=zeros(m-p,2*p);   RE2=zeros(m-p,2*p);
    else
        RE1=zeros(rex,p*p);   RE2=zeros(rex,p*p);
    end
    RF1=real(P);   RF2=imag(P);
    RM1=-real(KK);  RM2=-imag(KK);
    RN1=-real(B);   RN2=-imag(B);

    D11=kron(S2',eye1);
    DK=D11*K;  DS=D11*S;
    D=[DK,1i*DS];

    VAU=V'*A'*U;
    G11=kron(N,VAU);
    GK=G11*K;  GS=G11*S;
    G=[GK,1i*GS];

    VAUt=V'*A'*Ubot;
    H11=kron(N,VAUt);
    H=[H11,1i*H11];
    snx=0.5*p*(p-1); sny=2*p*(n-p); snx2=0.5*p*(p+1);
    if p==1
        SN1=zeros(p,2*(n-p));  SN2=zeros(p,2*(n-p));
    else
        SN1=zeros(snx,sny);  SN2=zeros(snx2,sny);
    end
    SM1=KK1*real(D);  SM2=KK2*imag(D);
    SE1=-KK1*real(G); SE2=-KK2*imag(G);
    SF1=-KK1*real(H); SF2=-KK2*imag(H);


    X11=kron(S2,eye4);
    X=[X11,1i*X11];

    VtAU=Vbot'*A'*U;
    Y11=kron(N,VtAU);
    YK=Y11*K;  YS=Y11*S;
    Y=[YK,1i*YS];

    VtAUt=Vbot'*A'*Ubot;
    Z11=kron(N,VtAUt);
    Z=[Z11,1i*Z11];
    tmx=p*(n-p); tmy=p*p;

    TN1=real(X);  TN2=imag(X);
    if p==1
        TM1=zeros(n-p,2*p);  TM2=zeros(n-p,2*p);
    else
        TM1=zeros(tmx,tmy);  TM2=zeros(tmx,tmy);
    end
    TE1=-real(Y); TE2=-imag(Y);
    TF1=-real(Z); TF2=-imag(Z);




    H=[QE1,QF1,QM1,QN1;
        QE2,QF2,QM2,QN2;
        RE1,RF1,RM1,RN1;
        RE2,RF2,RM2,RN2;
        SE1,SF1,SM1,SN1;
        SE2,SF2,SM2,SN2;
        TE1,TF1,TM1,TN1;
        TE2,TF2,TM2,TN2];


    % generate the right hand matrix g
    g1=real(U'*A*V*N);
    g2=imag(U'*A*V*N);
    g3=real(Ubot'*A*V*N);
    g4=imag(Ubot'*A*V*N);
    g5=real(V'*A'*U*N);
    g6=imag(V'*A'*U*N);
    g7=real(Vbot'*A'*U*N);
    g8=imag(Vbot'*A'*U*N);
    g=-[KK1*g1(:);
        KK2*g2(:);
        g3(:);
        g4(:);
        KK1*g5(:);
        KK2*g6(:);
        g7(:);
        g8(:)];


    % method 2 conjugate residual method
    [x,cell]=CR(H,g,itr);


    % construct Ex, Fx, Mx, Nx
    if p==1

        Ex1=0;
        Ex2=reshape(S*x(2), p,p);
        Ex=Ex1+1i*Ex2;
        Fx1=reshape(x(3:m+1), m-p,p);
        Fx2=reshape(x(m+2:2*m), m-p,p);
        Fx=Fx1+1i*Fx2;
        Mx1=0;
        Mx2=reshape(S*x(2*m+2), p,p);
        Mx=Mx1+1i*Mx2;
        Nx1=reshape(x(2*m+3:2*m+n+1), n-p,p);
        Nx2=reshape(x(2*m+n+2:2*m+2*n), n-p,p);
        Nx=Nx1+1i*Nx2;
    else
        pt=0.5*p*(p-1); mpt=m*p; nps=2*mpt-pp;

        Ex1=reshape(K*x(1:pt), p,p);
        Ex2=reshape(S*x(pt+1:pp), p,p);
        Ex=Ex1+1i*Ex2;

        Fx1=reshape(x(pp+1:mpt), m-p,p);
        Fx2=reshape(x(mpt+1:nps), m-p,p);
        Fx=Fx1+1i*Fx2;

        Mx1=reshape(K*x(nps+1:nps+pt), p,p);
        Mx2=reshape(S*x(nps+pt+1:2*mpt), p,p);
        Mx=Mx1+1i*Mx2;

        Nx1=reshape(x(2*mpt+1:2*mpt+n*p-pp), n-p,p);
        Nx2=reshape(x(2*mpt+n*p-pp+1:end), n-p,p);
        Nx=Nx1+1i*Nx2;
    end
    % construct xi and eta

    xi=U*Ex+Ubot*Fx;
    eta=V*Mx+Vbot*Nx;



    % genrate the next iterate V,P by cayley transform
    UG=U'*xi;
    VG=V'*eta;

    PZu=-xi+0.5*U*UG; PZv=-eta+0.5*V*VG;

    if invU
        WZu = PZu*U' - U*PZu';  Ru = WZu*U;
    else
        Uu=[PZu,U]; Wu=[U,-PZu];
        MU1=Wu'*U; MU2=Wu'*Uu;
    end
    WZv = PZv*V' - V*PZv';  Rv = WZv*V;
    if invU
        U = linsolve(eye(m) - 0.5*WZu, U + 0.5*Ru);
    else
        MMu=linsolve(eye2p-0.5*MU2,MU1);
        U=U+Uu*MMu;
    end
    V = linsolve(eye(n) - 0.5*WZv, V + 0.5*Rv);


    %-----终止判定标准------------
    [F,  Gu, Gv] = feval(fun, U, V, varargin{:});
    GU = Gu'*U; GV = Gv'*V;

    dtU = Gu - U*GU;     nrmGu  = norm(dtU, 'fro');
    dtV = Gv - V*GV;     nrmGv  = norm(dtV, 'fro');

    %---output data-------


    nrmG=sqrt(nrmGu^2+nrmGv^2);

    if (record >= 1)
        fprintf('%4d  %3.2e  %4.3e\n', ...
            itr, F, nrmG);
    end

    if  nrmG < 1e-6
        break;
    end


end
out.itr = itr;

function Kn = compute_Kn(n)

% 该程序给出了 当X 是反对称矩阵时，vec(X)=K_n vec_K(X)

Kn = [];
I = eye(n);
for k = 1:n-1
    Os = zeros(n*(k-1), n-k);
    Is = I(:, k+1:n);
    Ks = kron(eye(n-k), -I(:, k));
    Tmp = [Os; Is; Ks];
    Kn = [Kn Tmp];
end
Kn = Kn./sqrt(2);


function Sn = compute_Sn(n)
% 该程序给出了 当X 是对称矩阵时，vec(X)=S_n vec_S(X)

Sn = [];
I = eye(n);
for k = 1:n
    Os = zeros(n*(k-1), n-k+1);
    Is = I(:, k:n);
    Is(k,1) = sqrt(2);
    Ks = [zeros(n*(n-k), 1)  kron(eye(n-k), I(:,k))];
    Tmp = [Os; Is; Ks];
    Sn = [Sn Tmp];
end
Sn=Sn./sqrt(2);



function P = vecperm(m, n)
%VECPERM    Vec-permutation matrix.
%           VECPERM(M, N) is the vec-permutation matrix, an MN-by-MN
%           permutation matrix P with the property that if A is M-by-N then
%           vec(A) = P*vec(A').
%           If N is omitted, it defaults to M.

%   P is formed by taking every n'th row from EYE(M*N), starting with
%   the first and working down - see p. 277 of the reference.

%   Reference:
%   H. V. Henderson and S. R. Searle The vec-permutation matrix,
%   the vec operator and Kronecker products: A review Linear and
%   Multilinear Algebra, 9 (1981), pp. 271-288.

if nargin == 1, n = m; end

P = zeros(m*n);
I = eye(m*n);

k = 1;
for i=1:n
    for j=i:n:m*n
        P(k,:) = I(j,:);
        k = k+1;
    end
end