% Script for getting relevant and unique dimensionless groups 
clc; 
close all; 

% Choose parameter regime (1:Duffing, 2:2-DOF)
param_case = 2;
[q_l, q_u, caselabel, logcase, logflag] = set_param_case(param_case);

% First, do dimensional analysis.
%
% INPUTS
% m     = A(1), [1 0 0]
% c     = A(2), [1 0 -1]
% k1    = A(3), [1 0 -2]
% k3    = A(4), [1 -2 -2]
% D     = A(5), [2 0 -1]
%
% BY 
%
% M (kg), L (m), T (s)
%
% OUTPUTS
% Lyapunov, [0 0 -1]

D = [ 1 0 0; 1 0 0; 1 0 -1; 1 0 -1; 1 0 -2; 1 0 -2; 1 -1 -2; 2 0 -1; 2 0 -1]';
vq = [0 0 -1]';

% Set up
ndesign = 50;
fun = @(X) pressure_loss_colebrook(X);
[w, W, Q_ins, Q_outs] = alg1_rs_setup(D, vq, ndesign, fun, q_l, q_u, logflag);

% Run Algorithm 1
h = 1e-6;
nquad = 13;
[Z, e, U, C, delta, pi_outs] = alg1_rs( ndesign, w, W, Q_ins, Q_outs, h, nquad, q_l, q_u, logflag);
% save('D:\data-driven methods\系统参数灵敏度\二自由度非线性系统_Lyapunov\LyaData.dat','pi_outs','-ascii');

% Print out Z
Z

% With rotate variables
figure(1);
set(0,'defaultfigurecolor','w');
pp = get(gcf, 'PaperPosition');
pp(3)=5; pp(4)=5;
set(gcf, 'PaperPosition', pp);
scatter(delta(:, 2), delta(:, 3), 80, log(pi_outs), 'filled'); 
axis square; 
grid on;
set(gca, 'FontSize', 16);
% xlim([0.9*min(delta(:,1)) 1.1*max(delta(:,1))]);
xlabel('$\log({\pi}_2)$', 'interpreter', 'latex'); 
ylabel('$\log({\pi}_3)$', 'interpreter', 'latex');
title(sprintf('%s, %s', caselabel, logcase));
colorbar;
% print(sprintf('figs/ssp_alg1_%s', caselabel), '-depsc2');

% Express the new z variables in terms of the classical Reynolds number and
% roughness scale.
Ra = [0; 1; 1; 0; 1; -1; 0; 0; -1];
Rb = [1; 0; 0; 1; 1; -1; 0; 1; -1];
Rc = [-1; 1; 0; 0; 1; -1; 0; 1; -1];
Rd = [1; -1; 0; 0; 1; -1; 0; 0; 0];
Re = [-1; -1; 0; 1; 0; -1; 0; 2; -1];
Rf = [0; 1; 1; 0; 0; 0; 0; 0; -1];
W = [Ra  Rb  Rc  Rd  Re  Rf];

c = W \ Z(:,1);
fprintf('z_1 in terms of Ra Rb Rc Rd Re Rf: %6.4f, %6.4f, %6.4f, %6.4f, %6.4f, %6.4f\n', c(1), c(2), c(3), c(4), c(5), c(6) );
fprintf('Resid: %6.4e\n', norm(W*c - Z(:,1)) );
pai1 =Ra*c(1) +Rb*c(2) +Rc*c(3) +Rd*c(4) +Re*c(5) +Rf*c(6);
x1=sprintf('M%4.2f m%4.2f C%4.2f c%4.2f K%4.2f k%4.2f kk%4.2f D%4.2f d%4.2f',pai1(1),pai1(2),pai1(3),pai1(4),pai1(5),pai1(6),pai1(7),pai1(8),pai1(9))

c = W \ Z(:,2);
fprintf('z_2 in terms of Ra Rb Rc Rd Re Rf: %6.4f, %6.4f, %6.4f, %6.4f, %6.4f, %6.4f\n', c(1), c(2), c(3), c(4), c(5), c(6) );
fprintf('Resid: %6.4e\n', norm(W*c - Z(:,2)) );
pai2 =Ra*c(1) +Rb*c(2) +Rc*c(3) +Rd*c(4) +Re*c(5) +Rf*c(6);
x2=sprintf('M%4.2f m%4.2f C%4.2f c%4.2f K%4.2f k%4.2f kk%4.2f D%4.2f d%4.2f',pai2(1),pai2(2),pai2(3),pai2(4),pai2(5),pai2(6),pai2(7),pai2(8),pai2(9))

c = W \ Z(:,3);
fprintf('z_3 in terms of Ra Rb Rc Rd Re Rf: %6.4f, %6.4f, %6.4f, %6.4f, %6.4f, %6.4f\n', c(1), c(2), c(3), c(4), c(5), c(6) );
fprintf('Resid: %6.4e\n', norm(W*c - Z(:,3)) );
pai3 =Ra*c(1) +Rb*c(2) +Rc*c(3) +Rd*c(4) +Re*c(5) +Rf*c(6);
x3=sprintf('M%4.2f m%4.2f C%4.2f c%4.2f K%4.2f k%4.2f kk%4.2f D%4.2f d%4.2f',pai3(1),pai3(2),pai3(3),pai3(4),pai3(5),pai3(6),pai3(7),pai3(8),pai3(9))

c = W \ Z(:,4);
fprintf('z_4 in terms of Ra Rb Rc Rd Re Rf: %6.4f, %6.4f, %6.4f, %6.4f, %6.4f, %6.4f\n', c(1), c(2), c(3), c(4), c(5), c(6) );
fprintf('Resid: %6.4e\n', norm(W*c - Z(:,4)) );
pai4 =Ra*c(1) +Rb*c(2) +Rc*c(3) +Rd*c(4) +Re*c(5) +Rf*c(6);
x4=sprintf('M%4.2f m%4.2f C%4.2f c%4.2f K%4.2f k%4.2f kk%4.2f D%4.2f d%4.2f',pai4(1),pai4(2),pai4(3),pai4(4),pai4(5),pai4(6),pai4(7),pai4(8),pai4(9))

c = W \ Z(:,5);
fprintf('z_5 in terms of Ra Rb Rc Rd Re Rf: %6.4f, %6.4f, %6.4f, %6.4f, %6.4f, %6.4f\n', c(1), c(2), c(3), c(4), c(5), c(6) );
fprintf('Resid: %6.4e\n', norm(W*c - Z(:,5)) );
pai5 =Ra*c(1) +Rb*c(2) +Rc*c(3) +Rd*c(4) +Re*c(5) +Rf*c(6);
x5=sprintf('M%4.2f m%4.2f C%4.2f c%4.2f K%4.2f k%4.2f kk%4.2f D%4.2f d%4.2f',pai5(1),pai5(2),pai5(3),pai5(4),pai5(5),pai5(6),pai5(7),pai5(8),pai5(9))

c = W \ Z(:,6);
fprintf('z_6 in terms of Ra Rb Rc Rd Re Rf: %6.4f, %6.4f, %6.4f, %6.4f, %6.4f, %6.4f\n', c(1), c(2), c(3), c(4), c(5), c(6) );
fprintf('Resid: %6.4e\n', norm(W*c - Z(:,6)) );
pai6 =Ra*c(1) +Rb*c(2) +Rc*c(3) +Rd*c(4) +Re*c(5) +Rf*c(6);
x6=sprintf('M%4.2f m%4.2f C%4.2f c%4.2f K%4.2f k%4.2f kk%4.2f D%4.2f d%4.2f',pai6(1),pai6(2),pai6(3),pai6(4),pai6(5),pai6(6),pai6(7),pai6(8),pai6(9))


