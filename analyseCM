%----------------------------画图验证耦合矩阵的正确性
function [S11,S21]=analyseCM(n,Fstr,Fsto,M,CF,FBW,N_s)
% N,N_pole
% Fstr=Fd(1,1)
% Fsto=Fd(N_sample,1)
% M= coupling matrix
% N_s= N_sample
R=zeros(n+2,n+2);
R(1,1)=1;
R(n+2,n+2)=1;
W=eye(n+2,n+2);
W(1,1)=0;
W(n+2,n+2)=0;
for k=1:N_s
    wp(k)=((Fstr+(Fsto-Fstr)*k/N_s)/CF-CF/(Fstr+(Fsto-Fstr)*k/N_s))/FBW;
    Bd(k)=Fstr+(Fsto-Fstr)*k/N_s;
    A=R.*(-1i)+W.*wp(k)+M;
    AP=inv(A);
    S21(k,1)=-2*1i*AP(n+2,1);
    S11(k,1)=1+2*1i*AP(1,1);
end
