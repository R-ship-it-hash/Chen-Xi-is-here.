function  yLoss = J1_fun(x0)

D = [ 1 0 0; 1 0 0; 1 0 -1; 1 0 -1; 1 0 -2; 1 0 -2; 1 -1 -2; 2 0 -1; 2 0 -1]';
vq = [0 0 -1]';
% q_l = [-2 -2 -0.3 -0.3 -5 -5 0 -0.2 -0.2]; % lower parameter bound
% q_u = [4 4 0.4 0.4 10 10 5 0.4 0.4]; % upper parameter bound
h = 1e-20;

% Run Algorithm 1
I = [6 50 100 500 1000 5000 10000];
  for i=3:1:3
     ndesign = I(i);
     p1=x0(1);
     p2=x0(2);
     p3=x0(3);
     p4=x0(4);
     p5=x0(5);
     p6=x0(6);
     P = [p1;p2;p3;p4;p5;p6];

% Getting a basis for dimensional analysis
W = null(D);

% Make w consistent with the friction factor nondimensionalization
% w0 = D\vq;
w0 = lsqminnorm(D,vq);
w = w0 + W*P;   %dimensional
if norm(D*w-vq) > sqrt(eps), error('Bad normalization vector'); end
error0 = D*w-vq;

% Latin hypercube design
m = 9;
% X = lhsdesign(ndesign, m);
% Q = bsxfun(@plus,bsxfun(@times, X, (q_u-q_l)), q_l);  %将生成随机数调整到给定范围内
logflag = 0;
% 
% % Set bounds and draw quadrature points
% if logflag
%     % case of uniform density on logq
%     
%     logq_l = log(q_l); logq_u = log(q_u);
%     logQ = bsxfun(@plus,bsxfun(@times, 0.5*(X+1), (logq_u-logq_l)), logq_l);
%     Q_ins = logQ;
%     
% else
%     % case of uniform density on q
%     Q_ins = bsxfun(@plus,bsxfun(@times, 0.5*(X+1), (q_u-q_l)), q_l);
%     logQ = log(Q_ins);
%     
% end

A=xlsread('D:\data-driven methods\系统参数灵敏度\二自由度非线性系统_Lyapunov\parameter.xlsx');   %导入系统参数
logQ = log(A(1:100,1:9));
Q_ins = A(1:100,1:9);
Q_outs = A(1:100,12);
% if ndesign == 100
%    xlswrite('D:\data-driven methods\系统参数灵敏度\二自由度非线性系统_Lyapunov\parameter.xlsx',Q_ins);
% end

% Evaluate the pressure loss.
Q_outs = pressure_loss_colebrook(logQ);

% Set up
if logflag
    logQ = Q_ins;
else
    logQ = log(Q_ins);
end

% omega = quad_wts;
Nquad = length(Q_outs);

% Let's carefully follow the steps of Alg 2.

% (1) Compute evaluations of the dimensionless dependent variable
pi_outs = Q_outs.*exp(-logQ*w);

% (2) Compute evaluations of the logs of the dimensionless groups
gamma_ins = logQ*W;

% (3) Get perturbed points in logq space for finite differences
gamma_p = [gamma_ins(:,1)+h gamma_ins(:,2) gamma_ins(:,3) gamma_ins(:,4) gamma_ins(:,5) gamma_ins(:,6); gamma_ins(:,1) gamma_ins(:,2)+h gamma_ins(:,3) gamma_ins(:,4) gamma_ins(:,5) gamma_ins(:,6); gamma_ins(:,1) gamma_ins(:,2) gamma_ins(:,3)+h gamma_ins(:,4) gamma_ins(:,5) gamma_ins(:,6); gamma_ins(:,1) gamma_ins(:,2) gamma_ins(:,3) gamma_ins(:,4)+h gamma_ins(:,5) gamma_ins(:,6); gamma_ins(:,1) gamma_ins(:,2) gamma_ins(:,3) gamma_ins(:,4) gamma_ins(:,5)+h gamma_ins(:,6); gamma_ins(:,1) gamma_ins(:,2) gamma_ins(:,3) gamma_ins(:,4) gamma_ins(:,5) gamma_ins(:,6)+h];
logQ_p = (W' \ gamma_p')'; %200*9

% (4) For each finite difference point, run an experiment. By the way we
% chose the perturbed points, the first half of this vector contains
% perturbations in gamma_1 and the second half contains perturbations in
% gamma_2.
Q_outs_p = pressure_loss_colebrook(logQ_p);

% (5) Compute the corresponding dimensionless dependent variables
pi_outs_p = Q_outs_p.*exp(-logQ_p*w);

% (6) Compute finite differences
dg = [pi_outs_p(1:Nquad)-pi_outs pi_outs_p(Nquad+1:2*Nquad)-pi_outs pi_outs_p(2*Nquad+1:3*Nquad)-pi_outs pi_outs_p(3*Nquad+1:4*Nquad)-pi_outs pi_outs_p(4*Nquad+1:5*Nquad)-pi_outs pi_outs_p(5*Nquad+1:6*Nquad)-pi_outs] / h;
ndg = sqrt( sum(dg.^2, 2) );

% (7) Estimate active subspaces
% G = bsxfun(@times, dg, sqrt(omega)); %./ndg);
% C = G'*G;
G = dg ;
C = G'*G/ ndesign ;
% xlswrite('D:\data-driven methods\系统参数灵敏度\二自由度非线性系统_Lyapunov\C.xlsx',C);

% Compute eigenvalue decomposition and sort
[~, Sig, U] = svd(G, 'econ');
U;
e = diag(Sig).^2;
% fprintf('e1: %6.4e, e2: %6.4e, e3: %6.4e, e4: %6.4e, e5: %6.4e, e6: %6.4e\n', e(1), e(2), e(3), e(4), e(5), e(6));
yLoss = log10(e(2)) - log10(e(1));
% fprintf('FD rotates DA by %6.4f degrees\n', acos(U(1,1))*(180/pi));
Z = W*U;
  end
% fid=fopen(['D:\data-driven methods\系统参数灵敏度\二自由度非线性系统_Lyapunov\Jc','A.dat'],'a+');%写入文件路径
% for jj=1:length(yLoss)
% fprintf(fid,'%.4f\r\n',yLoss(jj));   %按列输出，若要按行输出：fprintf(fid,'%.4\t',A(jj)); 
% end
% fclose(fid);
end

