%%
% Script for running algorithms to get relevant and unique dimensionless
% groups. 
clc; 
clear all; 
tic
% Choose parameter regime (1:laminar, 2:turbulent, 3:high Re)
param_case = 2;
[q_l, q_u, caselabel, logcase, logflag] = set_param_case(param_case);

%%
% First, do dimensional analysis.
%
% INPUTS
% rho     = x(1), [1 -3 0]
% mu      = x(2), [1 -1 -1]
% D       = x(3), [0 1 0]
% epsilon = x(4), [0 1 0]
% U       = x(5), [0 1 -1]
%
% BY 
%
% M (kg), L (m), T (s)
%
% OUTPUTS
% dpdx, [1 -2 -2]

D = [ 1 0 0; 1 0 0; 1 0 -1; 1 0 -1; 1 0 -2; 1 0 -2; 1 -1 -2; 2 0 -1; 2 0 -1]';
% vq = [0 0 -1]';

% get the rational basis for comparison
W_null = null(D, 'r');

%%
% Set up
% nquad = 11;
% ndesign = 100;
% fun = @(X) pressure_loss_colebrook(X);
% [w, W, Q_ins, Q_outs] = alg2_fd_setup(ndesign, D, vq, fun, q_l, q_u, logflag);

%%
% Run Algorithm 2
% h = 1e-6;
% [Z, e, U, C, delta, pi_outs] = alg2_fd(ndesign, w, W, Q_ins, Q_outs, fun, h, logflag);

x0 = [0;0;0;0;0;0];
InF = 5;
A=[];
b=[];
Aeq=[];
beq=[];
lb=[-InF;-InF;-InF;-InF;-InF;-InF];
ub=[InF;InF;InF;InF;InF;InF];
a = fmincon('J1_fun',x0,A,b,Aeq,beq,lb,ub)
toc

%%

% % Let's make a 2-d summary plot 
% indz = randi(size(pi_outs,1), 1, 5000);
% 
% % With rotate variables
% figure(1);
% pp = get(gcf, 'PaperPosition');
% pp(3)=5; pp(4)=5;
% set(gcf, 'PaperPosition', pp);
% scatter(-delta(indz, 1), delta(indz, 2), 80, log(pi_outs(indz)), 'filled'); 
% axis square; grid on;
% set(gca, 'FontSize', 16);
% % xlabel('$\log(\hat{\pi}_1)$', 'interpreter', 'latex'); 
% % ylabel('$\log(\hat{\pi}_2)$', 'interpreter', 'latex');
% xlabel('$\hat{\pi}_1$', 'interpreter', 'latex'); 
% ylabel('$\hat{\pi}_2$', 'interpreter', 'latex');
% title(sprintf('%s, %s', caselabel, logcase));
% colorbar;
% 
% % With standard dimensionless groups (Re, epsilon/D)
% delta_da = log(Q_ins)*fliplr(W_null);
% 
% figure(2);
% pp = get(gcf, 'PaperPosition');
% pp(3)=5; pp(4)=5;
% set(gcf, 'PaperPosition', pp);
% scatter(delta_da(indz, 1), delta_da(indz, 2), 80, log(pi_outs(indz)), 'filled'); 
% axis square; grid on;
% set(gca, 'FontSize', 16);
% xlabel('$\pi_1$', 'interpreter', 'latex'); 
% ylabel('$\pi_2$', 'interpreter', 'latex');
% title(sprintf('%s, %s', caselabel, logcase));
% colorbar;


