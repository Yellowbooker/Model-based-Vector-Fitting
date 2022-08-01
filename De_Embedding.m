%% Microwave Filter Coupled-Resonator Circuit Modal Extraction By MVF
%% Version01----03/09/2020
%% Version02----03/27/2020
%% @Useless_Crap
%% First Step: Phase De-Embedding By VF
% Ref: Phase De-Embedding of Narrowband CoupledResonator Networks by
% Vector Fitting (Dr.Zhao Ping)

% % ======================= Clear Workspace  =======================
% 
% clear all
% clc
% 
% % ======================= Loop Parameters =======================
function [New_S11,New_S21,New_S22]=De_Embedding(F_band,S11,S21,S22,N,Nz,CF,BW)
X = 64;
err_ak_limit = 1e-3;
err_cp_limit = -30;	% in dB

% ======================= Phase Corrected  =======================

% ======================= Define Filter Parameter

% N = 6;
% Nz = 2;
% CF = 3498.6;
% BW = 200;
FBW = BW/CF;

% ======================= Define Loop Parameter

M = 2;
Line = 3;

% ======================= Bandpass to Lowpass
% 
% S_data = importdata('hfss6.s2p');
% F_band = S_data(:,1); % MHz
N_s = length(F_band); % number of sampling data
% S11(:,1) = S_data(:,2) + S_data(:,3)*1i; % S11
% S22(:,1) = S_data(:,8) + S_data(:,9)*1i; % S22
w_low = (F_band./CF - CF./F_band)/FBW;
s_low = 1i*w_low;

% ======================= Plot Sourse Data

% figure(10);
% plot(w_low,real(S11),'g--',w_low,imag(S11),'g','Linewidth',1.5);
% figure(11);
% plot(w_low,real(S22),'g--',w_low,imag(S22),'g','Linewidth',1.5);
% figure(12);
% PS22 = angle(S22);
% PS11 = angle(S11);
% plot(w_low,PS22,'y',w_low,PS11,'r--','LineWidth',1.5);

% =======================Constructing Systerm Matrices

M_S11 = zeros(N_s,N_s);
M_S22 = zeros(N_s,N_s);
A = zeros(N+M,N+M);
b = ones(N+M,1);
A1 = zeros(N_s,N+M+1);
A2 = zeros(N_s,N+M);

% ======================= Fit S11 and S22 Using VF =======================

% ======================= Initial Common Pole ak [-3,2.3]

for index = 1:N+M
    ak(index,1) = -3 + 5.3/(N + M - 1)*(index - 1);
    ak(index,1) = -0.01*abs(ak(index,1)) + ak(index,1)*1i;
end
tmp_ak = ak;
cp_plot = zeros(N+M, X);
err_ctl = zeros(1, X);

% ======================= Fit S11

for index = 1:X
    for index1 = 1:N_s
        for index2 = 1:N+M
            A2(index1,index2) = 1/(s_low(index1,1) - ak(index2,1));
        end
        M_S11(index1,index1) = S11(index1,1);
    end
    A1 = [A2,ones(N_s,1)];
    A = diag(ak,0);
    left = [A1   -1*M_S11*A2];
    right = S11;
%     Call = left\right;
    C_all = lsqminnorm(left,right);
    C = mat2cell(C_all,[N+M,1,N+M],1);
    cp = C{3,1};
    cp_plot(:,index)=10*log10(abs(cp));
    ak = eig(A-b*cp.');
    err_ctl(index) = sum(abs(ak - tmp_ak));
	if (err_ctl(index) < err_ak_limit) && (max(cp_plot(:,index)) < err_cp_limit)
		break
	end
	tmp_ak = ak;
end

%==============================	Disp loop time	======
disp('Loop time S11:')
disp(index)

% figure(1)
% subplot(1,2,1);
% plot([1:index],cp_plot([1:index]),'b','Linewidth',1.5);
% legend('cp', 'Location', 'NorthEast')
% title('Cp error in Log')
% grid on
% 
% subplot(1,2,2);
% plot([1:index],10*log10(err_ctl([1:index])),'r','Linewidth',1.5);
% legend('ak', 'Location', 'NorthEast')
% title('ak error in Log')
% grid on
% display(cp);
a11k = ak;
C11 = C{1,1};
d11 = C{2,1};

% ======================= Initial Common Pole ak [-3,2.3]

for index = 1:N+M
    ak(index,1) = -3 + 5.3/(N + M - 1)*(index - 1);
    ak(index,1) = -0.01*abs(ak(index,1)) + ak(index,1)*1i;
end
tmp_ak = ak;
cp_plot = zeros(N+M, X);
err_ctl = zeros(1, X);

% ======================= Fit S22

for index = 1:X
    for index1 = 1:N_s
        for index2 = 1:N+M
            A2(index1,index2) = 1/(s_low(index1,1) - ak(index2,1));
        end
        M_S22(index1,index1) = S22(index1,1);
    end
    A1 = [A2,ones(N_s,1)];
    A = diag(ak,0);
    left = [A1   -1*M_S22*A2];
    right = S22;
%     Call = left\right;
    C_all = lsqminnorm(left,right);
    C = mat2cell(C_all,[N+M,1,N+M],1);
    cp = C{3,1};
    ak = eig(A-b*cp.');
    cp_plot(:,index)=10*log10(abs(cp));
    ak = eig(A-b*cp.');
    err_ctl(index) = sum(abs(ak - tmp_ak));
	if (err_ctl(index) < err_ak_limit) && (max(cp_plot(:,index)) < err_cp_limit)
		break
	end
	tmp_ak = ak;
end
%==============================	Disp loop time	======
disp('Loop time S22:')
disp(index)

% figure(2)
% subplot(1,2,1);
% plot([1:index],cp_plot([1:index]),'b','Linewidth',1.5);
% legend('cp', 'Location', 'NorthEast')
% title('Cp error in Log')
% grid on
% 
% subplot(1,2,2);
% plot([1:index],10*log10(err_ctl([1:index])),'r','Linewidth',1.5);
% legend('ak', 'Location', 'NorthEast')
% title('ak error in Log')
% grid on
% display(cp);
a22k = ak;
C22 = C{1,1};
d22 = C{2,1};

% ======================= Caculate the Zeros of Sii

A11 = diag(a11k,0);
A22 = diag(a22k,0);
z11k = eig(A11-b*C11.'/d11);
z22k = eig(A22-b*C22.'/d22);

%---------------- Plot Zeros and Poles of S11 and S22

ReZ11 = real(z11k);
ImZ11 = imag(z11k);
ReP11 = real(a11k);
ImP11 = imag(a11k);
ReZ22 = real(z22k);
ImZ22 = imag(z22k);
ReP22 = real(a22k);
ImP22 = imag(a22k);
figure(3)
plot(ReZ11,ImZ11,'o',ReP11,ImP11,'x','LineWidth',1.5);
legend('S11 Zero', 'S11 Pole', 'NorthWest')
title('Zeroes and Poles of S11')
grid on
figure(4)
plot(ReZ22,ImZ22,'o',ReP22,ImP22,'x','LineWidth',1.5);
legend('S22 Zero', 'S22 Pole', 'NorthWest')
title('Zeroes and Poles of S22')
grid on

% ======================= Caculate S11 and S22 in Re and Im
TestNum = poly(z11k);
TestDen = poly(a11k);
for index = 1:N_s
    ExtrS11(index,1) = d11*polyval(TestNum,s_low(index,1))/polyval(TestDen,s_low(index,1));
end
TestNum = poly(z22k);
TestDen = poly(a22k);
for index = 1:N_s
    ExtrS22(index,1) = d22*polyval(TestNum,s_low(index,1))/polyval(TestDen,s_low(index,1));
end
% for index=1:N_s
%         for index1=1:N+M
%             trans1(index1,1) = s_low(index,1) - z11k(index1,1);
%             trans2(index1,1) = s_low(index,1) - a11k(index1,1);
%         end
%         ExtrS11(index,1) = d11*prod(trans1)/prod(trans2);
% end
% 
% for index=1:N_s
%         for index1=1:N+M
%             trans1(index1,1) = s_low(index,1) - z22k(index1,1);
%             trans2(index1,1) = s_low(index,1) - a22k(index1,1);
%         end
%         ExtrS22(index,1) = d22*prod(trans1)/prod(trans2);
% end

figure(5);
plot(w_low,real(S11),'y',w_low,imag(S11),'g',w_low,real(ExtrS11),'b--',w_low,imag(ExtrS11),'r--','Linewidth',1.5);
legend('Measured real(S11)', 'Measured imag(S11）','Fitted real(S11)', 'Fitted imag(S11）', 'NorthWest')
title('S11 in Real and Imag')
grid on
figure(6);
plot(w_low,real(S22),'y',w_low,imag(S22),'g',w_low,real(ExtrS22),'b--',w_low,imag(ExtrS22),'r--','Linewidth',1.5);
legend('Measured real(S22)', 'Measured imag(S22）','Fitted real(S22)', 'Fitted imag(S22）', 'NorthWest')
title('S22 in Real and Imag')
grid on

% ======================= Caculate Phase Factor
% ======================= S11
index2 = 0;
index3 = 0;
for index=1:N+M
        if abs(z11k(index,1))>=Line
                index2 = index2+1;
                alfa_z11k(index2,1) = z11k(index,1);
        end
        if abs(a11k(index,1))>=Line
                index3 = index3+1;
                alfa_a11k(index3,1) = a11k(index,1);
        end
end
TestNum = poly(alfa_z11k);
TestDen = poly(alfa_a11k);
for index = 1:N_s
    Extralfa11(index,1) = d11*polyval(TestNum,s_low(index,1))/polyval(TestDen,s_low(index,1));
end
% for index=1:N_s
%         for index1=1:length(alfa_a11k)
%             trans1(index1,1) = s_low(index,1) - alfa_z11k(index1,1);
%             trans2(index1,1) = s_low(index,1) - alfa_a11k(index1,1);
%         end
%         Extralfa11(index,1) = d11*prod(trans1)/prod(trans2);
% end
% figure(7)
% plot(w_low,abs(Extralfa11));
% figure(8)
% plot(w_low,angle(Extralfa11));
% figure(101);
% yyaxis left
% plot(w_low,abs(Extralfa11),'Linewidth',1.5);
% ylim([0,1.6]);
% yyaxis right
% plot(w_low,angle(Extralfa11),'Linewidth',1.5);
% legend('alfa in mag', 'alfa in phase', 'NorthWest')
% title('Alfa in S11')
% grid on

% ======================= S22
index2 = 0;
index3 = 0;
for index=1:N+M
        if abs(z22k(index,1))>=Line
                index2 = index2+1;
                alfa_z22k(index2,1) = z22k(index,1);
        end
        if abs(a22k(index,1))>=Line
                index3 = index3+1;
                alfa_a22k(index3,1) = a22k(index,1);
        end
end
TestNum = poly(alfa_z22k);
TestDen = poly(alfa_a22k);
for index = 1:N_s
    Extralfa22(index,1) = d22*polyval(TestNum,s_low(index,1))/polyval(TestDen,s_low(index,1));
end
% for index=1:N_s
%         for index1=1:length(alfa_a11k)
%             trans1(index1,1) = s_low(index,1) - alfa_z22k(index1,1);
%             trans2(index1,1) = s_low(index,1) - alfa_a22k(index1,1);
%         end
%         Extralfa22(index,1) = d22*prod(trans1)/prod(trans2);
% end

% figure(9)
% plot(w_low,abs(Extralfa22));
% figure(10)
% plot(w_low,angle(Extralfa22));

% ======================= Matlab R2016b Updata
% figure(102);
% yyaxis left
% plot(w_low,abs(Extralfa22),'Linewidth',1.5);
% ylim([0,1.6]);
% yyaxis right
% plot(w_low,angle(Extralfa22),'Linewidth',1.5);
% legend('alfa in mag', 'alfa in phase', 'NorthWest')
% title('Alfa in S22')
% grid on
% ======================= Unwrap Phase In Matlab
unwrap_Palfa11 = unwrap(angle(Extralfa11));
unwrap_Palfa22 = unwrap(angle(Extralfa22));
New_S11 = exp(-i.*angle(Extralfa11)).*S11;
New_S22 = exp(-i.*angle(Extralfa22)).*S22;
New_S21 = exp(-i.*(unwrap_Palfa11./2+unwrap_Palfa22./2)).*S21;
figure(9)
plot(w_low,angle(New_S11),'y',w_low,angle(New_S22),'r--','LineWidth',1.5);
% ======================= Polt Unwraped Phase
figure(102);
subplot(1,2,1);
yyaxis left
plot(w_low,abs(Extralfa11),'Linewidth',1.5);
ylim([0,1.6]);
yyaxis right
plot(w_low,unwrap_Palfa11,'Linewidth',1.5);
legend('alfa in mag', 'alfa in phase', 'NorthWest')
title('Alfa in S11')
grid on
subplot(1,2,2);
yyaxis left
plot(w_low,abs(Extralfa22),'Linewidth',1.5);
ylim([0,1.6]);
yyaxis right
plot(w_low,unwrap_Palfa22,'Linewidth',1.5);
legend('alfa in mag', 'alfa in phase', 'NorthWest')
title('Alfa in S22')
grid on
