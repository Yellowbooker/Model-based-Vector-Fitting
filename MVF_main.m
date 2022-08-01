%% Second Step: Modal-Based Vector Fitting (MVF) Technique
% Ref: Circuit Model Extraction of Parallel-Connected Dual-Passband
% Coupled-Resonator Filters (Dr.Zhao Ping)

% ======================= Clear Workspace =======================

clear all
clc

% ======================= Loop Parameters =======================

X = 64;
err_ak_limit = 1e-3;
err_cp_limit = -30;	% in dB

% ======================= Filter Parameters =======================

% Stop=3600;Start=3400;
% BW=Stop-Start;CF=sqrt(Start*Stop);FBW=BW/CF;
% N = 6;
% Nz = 4;
% CF = 3498.6;
% BW = 200;
% FBW = BW/CF;
N = 6;
Nz = 4;
CF = 1949.769217;
BW = 60;
FBW = BW/CF;
% N = 4;
% Nz = 0;
% CF = 29998.49996;
% BW = 9400;
% FBW = BW/CF;

% N = 10;
% Nz = 2;
% CF = 1048.81;
% BW = 100;
% FBW = BW/CF;

% N = 10;
% Nz = 4;
% CF = 2593.7666;
% BW = 162;
% FBW = BW/CF;
% N = 7;
% Nz = 2;
% CF = 1949.769217;
% BW = 60;
% FBW = BW/CF;

% ======================= Sampling data import =======================

S_data = importdata('hfss1.95.12.s2p');

F_band=S_data(:,1);%MHz
N_s=length(F_band);% number of sampling data
Origin_S11(:,1)=S_data(:,2)+S_data(:,3)*1i;%S11
Origin_S12(:,1)=S_data(:,4)+S_data(:,5)*1i;%S12
Origin_S21(:,1)=S_data(:,6)+S_data(:,7)*1i;%S21
Origin_S22(:,1)=S_data(:,8)+S_data(:,9)*1i;%S22
Origin_dBS11(:,1)=20*log10(abs(Origin_S11(:,1)));
Origin_dBS12(:,1)=20*log10(abs(Origin_S12(:,1)));
Origin_dBS21(:,1)=20*log10(abs(Origin_S21(:,1)));
Origin_dBS22(:,1)=20*log10(abs(Origin_S22(:,1)));
w_low=(F_band./CF-CF./F_band)/FBW;
s_low=1i*w_low;
% figure(1);
% plot(w_low,dBS11,'g--',w_low,dBS12,'g','Linewidth',1.5);
% legend('S11', 'S12', 'Location', 'NorthWest')
% title('Source S Parameter in Mag')
% grid on
% figure(2);
% PS22=angle(S22);
% PS11=angle(S11);
% plot(w_low,PS22,'b',w_low,PS11,'r','LineWidth',1.5);
% legend('S22', 'S11', 'Location', 'NorthWest')
% title('Measured S Parameter in Phase(rad/s)')
% grid on

[S11,S21,S22]=De_Embedding(F_band,Origin_S11,Origin_S21,Origin_S22,N,Nz,CF,BW);

S12=S21;

Z0=1;
for index=1:N_s
    [Y11(index,1),Y12(index,1),Y21(index,1),Y22(index,1)]=...
        StoY(Z0,S11(index,1),S12(index,1),S21(index,1),S22(index,1));
end
dBY11(:,1)=20*log10(abs(Y11(:,1)));
dBY12(:,1)=20*log10(abs(Y12(:,1)));
% % dBY21(:,1)=20*log10(abs(Y21(:,1)));
% % dBY22(:,1)=20*log10(abs(Y22(:,1)));
% figure(3);
% plot(w_low,dBY11,'r',w_low,dBY12,'b','Linewidth',1.5);
% legend('Y11', 'Y12', 'Location', 'NorthWest')
% title('Source Y Parameter in Mag')
% grid on


% ======================= Initial pole ak =======================

for index=1:N
    ak(index,1)=-1+2.1/(N-1)*(index-1);
    ak(index,1)=-0.01*abs(ak(index,1))+ak(index,1)*1i;
end

% ======================= Iterate Max X times =======================
M_Y11=zeros(N_s,N_s);
M_Y12=zeros(N_s,N_s);
M_Y22=zeros(N_s,N_s);
W12=zeros(N_s,N_s);
A=zeros(N,N);
Z1=zeros(N_s,Nz+1);
Z2=zeros(N_s,N+1);
b=ones(N,1);
A1=zeros(N_s,N+1);
A2=zeros(N_s,N);
A3=zeros(N_s,Nz+1);

tmp_ak = ak;
cp_plot = zeros(N, X);
err_ctl = zeros(1, X);

for index=1:X
    for index1=1:N_s
        for index2=1:N
            A2(index1,index2)=1/(s_low(index1,1)-ak(index2,1));
            trans(1,index2)=s_low(index1,1)-ak(index2,1);
        end
        for index3=1:Nz+1
            A3(index1,index3)=((s_low(index1,1))^(index3-1))/prod(trans);
        end
        E(index1,1)=1;
        M_Y11(index1,index1)=Y11(index1,1);
        M_Y12(index1,index1)=Y12(index1,1);
        M_Y22(index1,index1)=Y22(index1,1);
        W12(index1,index1)=1/sqrt(abs(Y12(index1,1)));
    end
    A1=[A2,E];
    for index1=1:N
        A(index1,index1)=ak(index1,1);
    end
%     [Q11,R11]=qr([A1,-M_Y11*A2]);
%     [Q12,R12]=qr([W12*A3,-W12*M_Y12*A2]);
%     [Q22,R22]=qr([A1,-M_Y22*A2]);
    left=[A1          Z1            Z2         -1*M_Y11*A2;
          Z2          W12*A3        Z2         -1*W12*M_Y12*A2;
          Z2          Z1            A1         -1*M_Y22*A2];
    right=[Y11;
           W12*Y12;
           Y22];
%     Call=left\right; % Caculate overdetermined equation in LS
    Call=lsqminnorm(left,right);
    C=mat2cell(Call,[N+1,Nz+1,N+1,N],1);
    cp=C{4,1};
    cp_plot(:,index)=10*log10(abs(cp));
    ak=eig(A-b*cp.');
%     display(ak);
%     display(cp);
	err_ctl(index) = sum(abs(ak - tmp_ak));
	if (err_ctl(index) < err_ak_limit) && (max(cp_plot(:,index)) < err_cp_limit)
		break
	end
	tmp_ak = ak;
end

%==============================	Disp loop time	======
disp('Loop time:')
disp(index)

figure(4)
subplot(1,2,1);
plot([1:index],cp_plot([1:index]),'b','Linewidth',1.5);
legend('cp', 'Location', 'NorthEast')
title('Cp error in Log')
grid on

subplot(1,2,2);
plot([1:index],10*log10(err_ctl([1:index])),'r','Linewidth',1.5);
legend('ak', 'Location', 'NorthEast')
title('ak error in Log')
grid on

C11=C{1,1};
C12=C{2,1};
C22=C{3,1};
% c11=A1\Y11;c22=A1\Y22;c12=(W12*A3)\(W12*Y12);

% ======================= Calculate ExtrY =======================
ExtrY11=zeros(N_s,1);
ExtrY21=zeros(N_s,1);
for index=1:N_s
    for index1=1:N
        ExtrY11(index,1)=ExtrY11(index,1)+(C11(index1,1)/(s_low(index,1)-ak(index1,1)));
        trans1(1,index1)=s_low(index,1)-ak(index1,1);
    end
    ExtrY11(index,1)=ExtrY11(index,1)+C11(N+1,1);
    ExtrY12(index,1)=polyval(flip(C12),s_low(index,1))/prod(trans1);
end
dBExtrY11(:,1)=20*log10(abs(ExtrY11(:,1)));
dBExtrY12(:,1)=20*log10(abs(ExtrY12(:,1)));
figure(5);
plot(w_low,dBY11,'y',w_low,dBY12,'g',w_low,dBExtrY11,'r--',w_low,dBExtrY12,'b--','Linewidth',1.5);
legend('Measured Y11', 'Measured Y12','Fitted Y11','Fitted Y12', 'Location', 'NorthWest')
title('Y Parameter in Mag')
grid on

% ======================= Find the R12 =======================
for index=1:N
    for index1=1:N-1
        if index1>=index
            flag=1;
        else
            flag=0;
        end
        trans2(1,index1)=ak(index,1)-ak(index1+flag,1);
    end
    R12(index,1)=polyval(flipud(C12),ak(index,1))/prod(trans2);
end       

% ======================= Find the Tzs =======================

Tz=roots(flipud(C12));
display(Tz);

% ======================= Y to TCM =======================

K11=C11(N+1,1);K22=C22(N+1,1);
M=zeros(N+2,N+2);
for index=1:N
    M(index+1,index+1)=1i*ak(index,1);
%     Q(index,1)=-1/imag(M(index+1,index+1))/FBW;
    if abs(C11(index,1))>=abs(C22(index,1))
        M(1,index+1)=sqrt(abs(C11(index,1)));
        M(N+2,index+1)=R12(index,1)/sqrt(abs(C11(index,1)));
    else if abs(C11(index,1))<abs(C22(index,1))
        M(N+2,index+1)=sqrt(abs(C22(index,1)));
        M(1,index+1)=R12(index,1)/sqrt(abs(C22(index,1)));
        end
    end
    M(index+1,1)=M(1,index+1);
    M(index+1,N+2)=M(N+2,index+1);
end
M(1,1)=K11/1i;
M(N+2,N+2)=K22/1i;
% display(M);
%% Third Step: Transform Matrix
% Ref: Microwave Filters for Communication Systems: Fundamentals, Design and Applications

% % ======================= folded coupling matrix =======================
% 
foldedCM=to_foldedCM(N,M);
% % disp('规范折叠拓扑：');
% display(real(foldedCM));
% 
% % ======================= wheel coupling matrix =======================
% 
wheelCM=to_wheelCM(N,foldedCM);
% % disp('规范轮形拓扑：');
% % display(wheelCM);
% 
% % ======================= CT coupling matrix =======================
% [CTCM,flag]=to_CTCM(N,Tz.',wheelCM);
% % disp('三元组（CT）拓扑：');
% % display(CTCM);
% 
% % ======================= CQ coupling matrix =======================
% 
% CQCM=to_CQCM(flag,N,CTCM,length(Tz));
% % disp('四元组（CQ）拓扑：');
% % display(CQCM);

for index=1:N
    Q(index,1)=-1/imag(foldedCM(index+1,index+1))/FBW;
end

figure(6)
bar(Q,'m');
title('Unloaded Quality Factor')
grid on

[ExtrS11,ExtrS21]=analyseCM(N,F_band(1,1),F_band(N_s,1),M,CF,FBW,N_s);
dBExtrS11(:,1)=20*log10(abs(ExtrS11(:,1)));
dBExtrS21(:,1)=20*log10(abs(ExtrS21(:,1)));
figure(7);
plot(F_band,Origin_dBS11,'y',F_band,Origin_dBS21,'g',F_band,dBExtrS11,'r--',F_band,dBExtrS21,'b--','Linewidth',1.5);
legend('Measured S11', 'Measured S21','Fitted S11','Fitted S21', 'Location', 'NorthWest')
title('S Parameter in Mag (T CM)')
xlim([F_band(1),F_band(N_s)]);
grid on

% %============================================================================================
% 
% [ExtrS11,ExtrS21]=analyseCM(N,F_band(1,1),F_band(N_s,1),foldedCM,CF,FBW,N_s);
% dBExtrS11(:,1)=20*log10(abs(ExtrS11(:,1)));
% dBExtrS21(:,1)=20*log10(abs(ExtrS21(:,1)));
% figure(8);
% plot(F_band,Origin_dBS11,'y',F_band,Origin_dBS21,'g',F_band,dBExtrS11,'r--',F_band,dBExtrS21,'b--','Linewidth',1.5);
% legend('Measured S11', 'Measured S21','Fitted S11','Fitted S21', 'Location', 'NorthWest')
% title('S Parameter in Mag (Folded CM)')
% grid on
% 
% %============================================================================================
% 
% [ExtrS11,ExtrS21]=analyseCM(N,F_band(1,1),F_band(N_s,1),CTCM,CF,FBW,N_s);
% dBExtrS11(:,1)=20*log10(abs(ExtrS11(:,1)));
% dBExtrS21(:,1)=20*log10(abs(ExtrS21(:,1)));
% figure(9);
% plot(F_band,Origin_dBS11,'y',F_band,Origin_dBS21,'g',F_band,dBExtrS11,'r--',F_band,dBExtrS21,'b--','Linewidth',1.5);
% legend('Measured S11', 'Measured S21','Fitted S11','Fitted S21', 'Location', 'NorthWest')
% title('S Parameter in Mag (CT CM)')
% grid on
% 
% %============================================================================================
% 
% [ExtrS11,ExtrS21]=analyseCM(N,F_band(1,1),F_band(N_s,1),CQCM,CF,FBW,N_s);
% dBExtrS11(:,1)=20*log10(abs(ExtrS11(:,1)));
% dBExtrS21(:,1)=20*log10(abs(ExtrS21(:,1)));
% figure(10);
% plot(F_band,Origin_dBS11,'y',F_band,Origin_dBS21,'g',F_band,dBExtrS11,'r--',F_band,dBExtrS21,'b--','Linewidth',1.5);
% legend('Measured S11', 'Measured S21','Fitted S11','Fitted S21', 'Location', 'NorthWest')
% title('S Parameter in Mag (CQ CM)')
% grid on
%% Output data in .txt

% fid=fopen(['G:\study\chen\MVF_V1_10032020\outputdata\','Extractive Matrix.txt'],'a');
% [r,c]=size(M);
% 
%  fprintf(fid,'\r\n');
%  fprintf(fid,'\r******************************************************\n');
%  fprintf(fid,'! Extractive Matrix Caculated by MVF\n');
%  fprintf(fid,'! Exported from MATLAB R2019b\n');
%  fprintf(fid,datestr(datetime));
%  fprintf(fid,'\r******************************************************\n');
%  fprintf(fid,'\rTCM:\n');
%  fprintf(fid,'\r\n');
%  for i=1:r
%    for j=1:c
%        if real(M(i,j)) < 0
%            fprintf(fid,'%f\t',real(M(i,j)));
%        else
%            fprintf(fid,'%f\t\t',real(M(i,j)));
%        end
%    end
%     fprintf(fid,'\r\n');
%  end
% fprintf(fid,'\r=================================================');
%  %========================================================
%  fprintf(fid,'\rfoldedCM:\n');
%  fprintf(fid,'\r\n');
%  for i=1:r
%    for j=1:c
%        if real(foldedCM(i,j)) < 0
%            fprintf(fid,'%f\t',real(foldedCM(i,j)));
%        else
%            fprintf(fid,'%f\t\t',real(foldedCM(i,j)));
%        end
%    end
%     fprintf(fid,'\r\n');
%  end
% fprintf(fid,'\r=================================================');
%  %========================================================
%  fprintf(fid,'\rWheelCM:\n');
%  fprintf(fid,'\r\n');
%  for i=1:r
%    for j=1:c
%        if real(wheelCM(i,j)) < 0
%            fprintf(fid,'%f\t',real(wheelCM(i,j)));
%        else
%            fprintf(fid,'%f\t\t',real(wheelCM(i,j)));
%        end
%    end
%     fprintf(fid,'\r\n');
%  end
% fprintf(fid,'\r=================================================');
%  %========================================================
%  fprintf(fid,'\rCTCM:\n');
%  fprintf(fid,'\r\n');
%  for i=1:r
%    for j=1:c
%        if real(CTCM(i,j)) < 0
%            fprintf(fid,'%f\t',real(CTCM(i,j)));
%        else
%            fprintf(fid,'%f\t\t',real(CTCM(i,j)));
%        end
%    end
%     fprintf(fid,'\r\n');
%  end
% fprintf(fid,'\r=================================================');
%  %========================================================
%  fprintf(fid,'\rCQCM:\n');
%  fprintf(fid,'\r\n');
%  for i=1:r
%    for j=1:c
%        if real(CQCM(i,j)) < 0
%            fprintf(fid,'%f\t',real(CQCM(i,j)));
%        else
%            fprintf(fid,'%f\t\t',real(CQCM(i,j)));
%        end
%    end
%     fprintf(fid,'\r\n');
%  end
% fprintf(fid,'\r=================================================');
