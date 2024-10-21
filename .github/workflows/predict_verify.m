clc
clear all
% A=xlsread('D:\data-driven methods\系统参数灵敏度\二自由度非线性系统_Lyapunov\difference_step.xlsx');   
% B=A(21:27,1:7);
% C=B';

Pi=load('LyaPi.dat');
Data=load('LyaData.dat');
x=(1:1:length(Pi))';
plot(x,Pi);
hold on;
plot(x,Data);
legend('Pi','Data');


