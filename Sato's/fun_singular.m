function [F,Gu,Gv]=fun_singular(U,V,A,N)
%% ------- ������������ֵ�ֽ�����--------------------------------
%  min  F(U,V)= -Re(tr(U'AVN))  (U,V) \in St(m,p,C)*St(n,p,C)
%                                U'*U = I_p, where U \in C^{m,p}
%                                V'*V = I_p, where V \in C^{n,p}
%% -------------------------------------------------------------------------
UA=U'*A*V*N;
F=-trace(UA); % objective function
Gu=-A*V*N;       % the derivation of  f(x) with respect to V    ŷʽ�ݶ�
Gv=-A'*U*N;      % the derivation of  f(x) with respect to P    ŷʽ�ݶ�
