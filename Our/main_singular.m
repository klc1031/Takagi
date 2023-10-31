clc
clear
close all

%% For Matlab solver for Riemannian optimization on the Stiefel manifold
%  min -1/2 tr(Jp U'*Jn*A*U*N)
%  s.t.   U'*U = I2p


opts.record = 0;
opts.mxitr  = 3000;
opts.xtol = 1e-6;
opts.gtol = 1e-6;
opts.ftol = 1e-12;
opts.mm = 2;  % backward integer for nonmonotone line search
opts.nt = 5;
opts.alpha = 1.4;
opts.lambda = 0.2;  % factor for decreasing the step size in the backtracking line search
opts.delta = 1e-4;  % parameters for control the linear approximation in line search

W=[4,20,50,100,400,100,300,500, 300];
P=[2,  2,2,  2,  2,  3,  3,  3,  5];

for i= 1
    n=W(i);p=P(i);

    %% Example 5.1
    c=[1+2i,2-4i,2-1i,-1+3i];
    r=[-1+3i,3-1i,1-2i,4-2i];
    A1 = hankel(c,r);

    %% Example 5.2
    %     A1 = randn(n,n) + 1i*randn(n,n);
    %     c1=randn(1,n-1) + 1i*randn(1,n-1);
    %     c = [c1, 1];
    %     r1=randn(1,n-1) + 1i*randn(1,n-1);
    %     r = [1, r1];
    %     A1 = hankel(c,r);

    %% real form
    A = [real(A1),imag(A1);-imag(A1),real(A1)];
    N1=[p:-1:1];
    N2=[N1,N1];
    N=diag(N2);
    U1 = randn(n,p)+1i*randn(n,p);
    U2 = orth(U1);
    U0 = [real(U2),imag(U2);-imag(U2),real(U2)];
    B1=eye(n);
    C1=zeros(n);
    J1=[C1,B1;B1,C1];
    B2=eye(p);
    C2=zeros(p);
    J2=[C2,B2;B2,C2];

    figure(i)
    hold on
    %% 最速下降法
    tic;
    [out]= wen_Stiefel_singular(U0, @fun_singular, opts, A,N,J1,J2);
    Time_wen = toc;
    Witr=out.itr;
    Wfval=out.fval;
    WGrad=out.gradfu;
    TTwen=out.TT;
    U_1=out.U;

    %%    基于QR分解的共轭梯度法
    tic;
    [out]= zhu_StiQR_singular(U0, @fun_singular, opts, A,N,J1,J2);
    Citr=out.itr;
    Time_zhu = toc;
    zfval=out.fval;
    zGrad=out.Gradu;
    TTCG=out.TT;
    U_2=out.U;

    %% Newton 法
    tic;
    [U_new,TT_newton,itr_new,Fu_new,Grad_new]=Newton(U0, @fun_singular,opts,A,N,J1,J2);
    Time_new=toc;
    Nitr1=itr_new;

    %% 牛顿法（混合）
    tic;
    [out]= zhu_StiQR_singular2(U0, @fun_singular, opts, A, N,J1,J2);
    itr_1=out.itr;
    U_h1=out.U;
    TT_1=out.TT;
    [U_hy,G_hy,itr2,F_hy,G_h2]=Newton(U_h1, @fun_singular,opts,A,N,J1,J2);
    time_hybrid=toc;
    itr_2=itr2;
    Nitr=itr_1+itr_2;
    TT_2=G_hy;
    TT_hynewton=[TT_1,TT_2];

    %% print numerical results
    fprintf('%6s\t %4s\t %4s\t %4s\t %4s\t %4s\t\n', ...
        'Method','    [n,p]',' CPU','IT','  Fvalue','Grad');
    fprintf('%4s\t %2d %2d\t  %1.3f\t %1.3f\t   %4.2e\t %4.2e\n', ...
        'alg4.1', n,p, Time_wen,Witr, Wfval,WGrad);
    fprintf('%4s\t %2d %2d\t %1.3f\t %1.3f\t   %4.2e\t %4.2e\n', ...
        'alg4.2', n,p, Time_zhu,Citr, zfval,zGrad);
    fprintf('%4s\t %2d %2d\t %1.3f\t %1.3f\t   %4.2e\t %4.2e\n', ...
        'alg4.3', n,p, Time_new,Nitr1,Fu_new,Grad_new);
    fprintf('%4s\t %2d %2d\t %1.3f\t %1.3f\t   %4.2e\t %4.2e\n', ...
        'alg4.4', n,p, time_hybrid, Nitr,F_hy,G_h2);
    fprintf('---------------------------------------------------------------------------------------- \n');

    %% figure
    if Witr<300
        h_Wzhu=plot(log10(TTwen),  'm+-.','Markersize',6,'LineWidth',1);
    else
        h_Wzhu=plot(log10(TTwen(1:300)),  'm+-.','Markersize',6,'LineWidth',1);
    end

    if Citr<300
        h_zhuQR=plot(log10(TTCG),  'b^-.','Markersize',4,'LineWidth',0.6);
    else
        h_zhuQR=plot(log10(TTCG(1:300)),  'b^-.','Markersize',4,'LineWidth',0.6);
    end

    if Nitr1<300

        h_new=plot(log10(TT_newton),  'ro:','Markersize',6,'LineWidth',1);
    else
        h_new=plot(log10(TT_newton(1:300)),  'ro:','Markersize',6,'LineWidth',1);
    end

    if Nitr<300
        h_hynewton=plot(log10(TT_hynewton),  'cd-.','Markersize',6,'LineWidth',1);
    else
        h_hynewton=plot(log10(TT_hynewton(1:300)),  'cd-.','Markersize',6,'LineWidth',1);
    end

    xlabel('Iteration','FontSize',15);
    ylabel('$\|gradF(U)\|$','Interpreter','latex','fontsize',10)
    yticks(-12:2:4);
    yticklabels({'10^{-12}','10^{-10}','10^{-8}','10^{-6}','10^{-4}','10^{-2}','10^{0}','10^{-2}','10^{4}',});
    legend('Alg4.1','Alg4.2','Alg4.3','Alg4.4')
    set(gca,'LineWidth',1)
    set(gca,'FontSize',15)
    title(['n=',num2str(n),',p=',num2str(p)],'fontsize',14)
    hold off
end