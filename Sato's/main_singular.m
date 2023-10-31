clc
clear
close all

m=4; n=4; p=2;
opts.record = 0;
opts.mxitr  = 3000;
opts.xtol = 1e-6;
opts.gtol = 1e-6;

opts.ftol = 1e-12;
opts.mm = 2;  % backward integer for nonmonotone line search
opts.nt = 5;
opts.alpha = 1e-3;
opts.lambda = 0.2;  % factor for decreasing the step size in the backtracking line search
opts.delta = 1e-4;  % parameters for control the linear approximation in line search


%% 问题模型为 复奇异值分解问题
%  A=U*sigma*V^H     A \in K^(m,n)   U \in U(m),   V\in U(n)

%  问题转化为最小迹函数问题
%  min  F(U,V)= -Re(tr(U'AVN))  (U,V) \in St(m,p,C)*St(n,p,C)
%                                U'*U = I_p, where U \in C^{m,p}
%                                V'*V = I_p, where V \in C^{n,p}



c=[1+2i,2-4i,2-1i,-1+3i];
r=[-1+3i,3-1i,1-2i,4-2i];
A = hankel(c,r);
N1=[2:-1:1];
N=diag(N1);


%% 初始未知量
U0 = randn(m,p)+1i*randn(m,p);    U0 = orth(U0);
V0 = randn(n,p)+1i*randn(n,p);    V0 = orth(V0);
rng("default") %固定初值的选取，可删除

%% 最速下降法
tic;
[out]= Alg1(U0, V0, @fun_singular, opts, A,N);
Witr=out.itr;
Time_wen = toc;


%% 共轭梯度法
tic;
[ out]= Alg2(U0, V0, @fun_singular, opts, A,N);
Time_zhu = toc;
zitr=out.itr;

%% 牛顿法
[U,V,out]=Alg3(U0, V0, @fun_singular, opts, A,N);
time_new=toc;
Newitr=out.itr;

%% 牛顿法（混合）
tic;
opts.tt = 1;
opts.gtol = 0.1;
[out]= Alg2(U0, V0, @fun_singular, opts, A,N);
U0=out.U;
V0=out.V;
itr_1=out.itr;
[U,V,out]=Alg3(U0, V0, @fun_singular, opts, A,N);
time_hybrid=toc;
itr_2=out.itr;
Nitr=itr_1+itr_2;



fprintf('%6s\t  %4s\t %4s\t \n', ...
    'Method',' CPU','IT');
fprintf('%4s\t  %1.3f\t %3.1f\n', ...
    'Alg1',Time_wen,Witr);
fprintf('%4s\t  %1.3f\t %3.1f\n', ...
    'Alg2', Time_zhu,zitr);
fprintf('%4s\t  %1.3f\t %3.1f\n', ...
    'Alg3', time_new, Newitr); %%%%牛顿法
fprintf('%4s\t  %1.3f\t %3.1f\n', ...
    'Alg4 ',  time_hybrid, Nitr);