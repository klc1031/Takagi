function [M]=xishu(n,p,U,A,N,Jn,Jp)
% The coefficient matrix M for solving Newton equation

%% -------------------------------------------------------------------------
S=U'*Jn*A*U*Jp*N;
S1=0.5*(S+S');
E1=eye(n);
M1=kron(S1,E1);
S2=Jp*N;
E2=0.5*U*U'*Jn*A-Jn*A;
M2=kron(S2,E2);
E3=U*U';
M3=-0.5*kron(S1,E3);
S4=eye(p);
E4=U*S1*U';
M4=0.5*kron(S4,E4);

M50=zeros(n*p,n*p);
S5=U'*Jn*A;
E5=U*Jp*N;
M51=0.5*kron(S5,E5);
for i=1:n
    for j=1:p
        S5_1=zeros(n,p);
        S5_1(i,j)=1;
        S5_2=S5_1;

        M52=kron(S5_2,S5_2');

        M50=M50+M52;
    end
end
M5=M51*M50;
M=M1+M2+M3+M4+M5;
