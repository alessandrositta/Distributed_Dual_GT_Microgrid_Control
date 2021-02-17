%
% Distributed Control System
% Course Project 2020 
% Distributed Dual Gradient Tracking for Microgrid Control
%
% Task 1
% 
% Bologna, 24/01/2020
%

%% Initialization
clear all;
clc; close all;
rng(1)

NN = 5;

p = 0.1;
I_NN = eye(NN,NN);
notI = ~I_NN;

%% Graph creation

while 1
  Adj = binornd(1,p,NN,NN); % Randomly generated matrix
  
  Adj = Adj.*notI;    % Set diagonal to 0 (remove self-loops)
  Adj = or(Adj,Adj'); % Symmetrize

  % Test connection
  test = (I_NN+Adj)^NN;
  if ~any(any(~test))
    fprintf('the graph is connected\n');
    break
  else
    fprintf('the graph is not connected\n');
  end
end

% Doubly stochastic

DEGREE = sum(Adj);
WW = zeros(NN,NN);

for ii = 1:NN
  N_ii = find(Adj(:,ii) == 1)';
  for jj = N_ii
    WW(ii,jj) = 1/(1 + max(DEGREE(ii),DEGREE(jj) ));
  end
  WW(ii,ii) = 1 - sum(WW(ii,:));
end

%% Quadratic cost parameters initialization

% q_i(x) = - (x^T Q_i x + R_i^T x)
% grad_q_i(x) = - (2*(Q_i+Q_i^T) x + R_i)
% we change the problem to a minimazion one considering -q(x) instead of
% q(x)

SS = 3; % x in R^S

QQ = zeros(SS,SS,NN);
RR = zeros(SS,NN);
a = 1;   b = 10;
c = -5;  d = 5;
for ii = 1:NN
  diag_matrix = diag((b-a)*rand(SS,1) + a); % entries unif. in [a,b]
  T = orth(rand(SS,SS)); % change of coordinates
  QQ(:,:,ii) = T*diag_matrix*T';
  QQ(:,:,ii) = (QQ(:,:,ii)'+QQ(:,:,ii))/2; % Symmetrize
  
  RR(:,ii) = (d-c)*rand(SS,1) + c;  % entries unif. in [c,d]
end

%% Centralized solution

QQ_centr = sum(QQ,3); %summing over the last column the elements of QQ
RR_centr = sum(RR,2);

% min_x q = - \sum q_i =  x^T (sum_i Q_i) x + ( sum_i R_i)^T*x
[xopt,qopt] = quadprog(2*QQ_centr,RR_centr);

%% Distributed setup

MAXITERS = 1e4;

XX = zeros(SS,NN,MAXITERS);
YY = zeros(SS,NN,MAXITERS);

values        = zeros(MAXITERS,1);
consensus_err = zeros(MAXITERS,1);
tracking_err  = zeros(MAXITERS,1);

x_init  = rand(SS,NN);
XX(:,:,1) = x_init;

% initialize the tracker YY
for ii=1:NN
  [~,YY(:,ii,1)] = quadratic_function(XX(:,ii,1),QQ(:,:,ii),RR(:,ii));
end
% initialization of a constant stepsize
gamma = 1*1e-3;

for tt = 1:MAXITERS-1
  fprintf('iteration %d\n',tt);

  % Solution estimate update
  for ii=1:NN
    N_ii = find(Adj(:,ii) == 1)';

    % x_i^t+1 = sum_{j \in N_i U { i } } w_ij x_j^t - gamma* y_i^t
    XX(:,ii,tt+1) = WW(ii,ii)*XX(:,ii,tt);

    for jj = N_ii
      XX(:,ii,tt+1) = XX(:,ii,tt+1) + WW(ii,jj)*XX(:,jj,tt);
    end
    
    XX(:,ii,tt+1) = XX(:,ii,tt+1) - gamma*YY(:,ii,tt);
  end
  
  % Tracker update
  for ii=1:NN
    N_ii = find(Adj(:,ii) == 1)';

    YY(:,ii,tt+1) = WW(ii,ii)*YY(:,ii,tt);

    for jj = N_ii
      YY(:,ii,tt+1) = YY(:,ii,tt+1) + WW(ii,jj)*YY(:,jj,tt);
    end

    [~,grad_f_ii_new] = quadratic_function(XX(:,ii,tt+1),QQ(:,:,ii),RR(:,ii));
    [~,grad_f_ii_old] = quadratic_function(XX(:,ii,tt),  QQ(:,:,ii),RR(:,ii));

    YY(:,ii,tt+1) = YY(:,ii,tt+1) + (grad_f_ii_new-grad_f_ii_old);
  end
  
  % Performance check
  XX_avg = mean(XX(:,:,tt),2);
  YY_avg = mean(YY(:,:,tt),2);

  for ii=1:NN
    [q_ii,~] = quadratic_function(XX(:,ii,tt),QQ(:,:,ii),RR(:,ii));
    values(tt) = values(tt) + q_ii;
    consensus_err(tt) = consensus_err(tt) + norm(XX(:,ii,tt) - XX_avg);
    tracking_err(tt)  = tracking_err(tt)  + norm(YY(:,ii,tt) - YY_avg);
  end
end

% Last value [for plot]
tt = MAXITERS;
XX_avg = mean(XX(:,:,tt),2);
YY_avg = mean(YY(:,:,tt),2);

for ii=1:NN
  q_ii = quadratic_function(XX(:,ii,tt),QQ(:,:,ii),RR(:,ii));
  values(tt) = values(tt) + q_ii;
  consensus_err(tt) = consensus_err(tt) + norm(XX(:,ii,tt) - XX_avg);
  tracking_err(tt)  = tracking_err(tt)  + norm(YY(:,ii,tt) - YY_avg);
end

%% Plot of the results

figure
  semilogy(1:MAXITERS,abs(values(1:MAXITERS)-qopt), 'LineWidth',2);
  hold on, grid on
  xlabel('t')
  ylabel('cost error')

figure
  semilogy(1:MAXITERS,consensus_err(1:MAXITERS), 'LineWidth',2);
  hold on, grid on
  xlabel('t')
  ylabel('consensus error')

figure
  semilogy(1:MAXITERS,tracking_err(1:MAXITERS), 'LineWidth',2);
  hold on, grid on
  xlabel('t')
  ylabel('tracking error')
  
  
figure
  plot(1:MAXITERS, squeeze(XX(1,:,1:MAXITERS)),'LineWidth',2);
  hold on, grid on
  plot(1:MAXITERS, xopt(1)*ones(MAXITERS,1),'--','LineWidth',2);
  hold on, grid on
  plot(1:MAXITERS, squeeze(XX(2,:,1:MAXITERS)),'LineWidth',2);
  hold on, grid on
  plot(1:MAXITERS, xopt(2)*ones(MAXITERS,1),'--','LineWidth',2);
  hold on, grid on
  plot(1:MAXITERS, squeeze(XX(3,:,1:MAXITERS)),'LineWidth',2);
  hold on, grid on
  plot(1:MAXITERS, xopt(3)*ones(MAXITERS,1),'--','LineWidth',2);
  xlabel('t')
  ylabel('x_i^t, x^*') 
  
  
 