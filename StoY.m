function [Y11,Y21,Y12,Y22]=StoY(Z0,S11,S12,S21,S22)
%针对滤波器的参数转化
Y11=((1-S11)*(Z0+S22*Z0)+S21*S12*Z0)/((Z0+S11*Z0)*(Z0+S22*Z0)-S12*S21*Z0*Z0);
Y12=-2*S12*Z0/((Z0+S11*Z0)*(Z0+S22*Z0)-S12*S21*Z0*Z0);
Y21=Y12;
Y22=((1-S22)*(Z0+S11*Z0)+S21*S12*Z0)/((Z0+S11*Z0)*(Z0+S22*Z0)-S12*S21*Z0*Z0);