function [F,Gu]=fun_singular(U,A,N,J1,J2)
%% -------------------------------------------------------------------------
UA=J2*U'*J1*A*U*N;
F=-0.5*trace(UA); % objective function
Gu=-0.5*(A'*J1+J1*A)*U*J2*N; % the derivation of  F(U) with respect to U   

