function [F,Gu,Gv]=fun_singular(U,V,A,N)
%% ------- 复数域上奇异值分解问题--------------------------------
%  min  F(U,V)= -Re(tr(U'AVN))  (U,V) \in St(m,p,C)*St(n,p,C)
%                                U'*U = I_p, where U \in C^{m,p}
%                                V'*V = I_p, where V \in C^{n,p}
%% -------------------------------------------------------------------------
UA=U'*A*V*N;
F=-trace(UA); % objective function
Gu=-A*V*N;       % the derivation of  f(x) with respect to V    欧式梯度
Gv=-A'*U*N;      % the derivation of  f(x) with respect to P    欧式梯度
