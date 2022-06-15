clc
clear

R=[];
A=[];
% B=[];
% C=[];
% D=[];
% E=[];



load('D:\ZN_Matlab_Calculate\ROC_AUC_Calculate\true_result')
R=true_result;
load('D:\ZN_Matlab_Calculate\ROC_AUC_Calculate\predict1_result')
A=predict1_result;
% load('D:\MATLAB--Calculate\Feature_Cauculate\predict2_result')
% B=predict2_result;
% load('D:\MATLAB--Calculate\Feature_Cauculate\predict3_result')
% C=predict3_result;
% load('D:\MATLAB--Calculate\Feature_Cauculate\predict4_result')
% D=predict4_result;
% load('D:\MATLAB--Calculate\Feature_Cauculate\predict5_result')
% E=predict5_result;

Q1=value(A,R);
% Q2=value(B,R);
% Q3=value(C,R);
% Q4=value(D,R);
% Q5=value(E,R);
% Q=[Q1;Q2;Q3;Q4;Q5];
% S1=[1;2;3;4;5];
% S2={'value','TP','TN','FP','FN','acc','pre','sen','spe','mcc','f'};
% S3=[S1,Q1];
% S4=num2cell(S3);


S2={'TP','TN','FP','FN','acc','pre','sen','spe','mcc','f'};
S4=num2cell(Q1);
S=[S2;S4];
xlswrite('predict_data2.xls',S);