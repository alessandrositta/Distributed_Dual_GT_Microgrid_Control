% Distributed Control System
% Course Project 2020 
% Distributed Dual Gradient Tracking for Microgrid Control
%
% Task 2
% Version 10
% 
% Bologna, 10/02/2020
%
clear all;
clc; close all;
rng(1)

%% General initializations

GG = 4; %number of GENERATORS in the Network
SS = 3; %number of STORAGES in the Network
CC = 2; %number of CONTROLLABLE_LOADs in the Network
NN = GG+SS+CC+1; %total number of nodes in the Network
TT = 12; % p_i in R^(T+1)

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
%% Definition of the Network

% Matrices for Centralized Solution
AA_tot = [];
bb_tot = [];
QQ_centr = [];
RR_centr = [];

% Matrices for the Cost Functions
QQ_g =[];
RR_g =[];
QQ_s =[];
RR_s =[];
QQ_c =[];
RR_c =[];

% Matrices for the Local Constraints
AA_g =[];
bb_g =[];
AA_s =[];
bb_s =[];
AA_c =[];
bb_c =[];

%% GENERATORS parameters initialization

% [AA_tot,bb_tot,QQ_centr,RR_centr] = generator_param(GG,alpha1,alpha2,ub_1_g,lb_2_g,ub_3_g,lb_4_g);
[AA_tot,bb_tot,QQ_centr,RR_centr,AA_g,bb_g,QQ_g,RR_g] = generator_param(TT,GG,1,1,5,-1,5,-1,AA_tot,bb_tot,QQ_centr,RR_centr);

%% STORAGEs parameters initialization

% [AA_tot,bb_tot,QQ_centr,RR_centr] = storage_param(SS,epsilon,ub_1_s,lb_2_s,ub_3_s,lb_4_s);
[AA_tot,bb_tot,QQ_centr,RR_centr,AA_s,bb_s,QQ_s,RR_s] = storage_param(TT,SS,1,3,-0.5,4,0,AA_tot,bb_tot,QQ_centr,RR_centr);

%% CONTROLLABLE_LOADs parameters initialization

% [AA_tot,bb_tot,QQ_centr,RR_centr,AA_c,bb_c,QQ_c,RR_c] = conl_param(TT,CC,epsilon,beta,P_des,P_bound,AA_tot_in,bb_tot_in,QQ_centr_in,RR_centr_in)
time = [0:TT-1];
P_des = [];
for ii=1:CC
    P_des = [P_des,(sin(time)+1*rand)'];
end

[AA_tot,bb_tot,QQ_centr,RR_centr,AA_c,bb_c,QQ_c,RR_c] = conl_param(TT,CC,2,3,P_des,5,AA_tot,bb_tot,QQ_centr,RR_centr);

%% POWER TRADING COST parameter initialization

%[AA_tot,bb_tot,QQ_centr,RR_centr,AA_t,bb_t,QQ_t,RR_t] = trading_param(TT,epsilon,c1,c2,E,AA_tot_in,bb_tot_in,QQ_centr_in,RR_centr_in)
[AA_tot,bb_tot,QQ_centr,RR_centr,AA_t,bb_t,QQ_t,RR_t] = trading_param(TT,1,2,3,5,AA_tot,bb_tot,QQ_centr,RR_centr);

%% Centralized solution
options_qp = optimoptions('quadprog','Display','none');

% min_p f =  \sum f_i(p_i) = p_i^T Q_i x + R_i^T p_i)
% subject to [1...1]*[p_1 ... p_N]^T = 0 (globabl constraints)
% p_i belonging to X_i (local constraints)

%definition of the EQUALITY constraint (in our case: coupling
%constraints)

% Construction of AA_eq for the GENERATORS
aux_r_g = zeros(1,2*GG*TT);
aux_r_g(1,1)=1;
for ii=1:(GG-1)
    aux_r_g(1,(ii*(2*TT))+1)=1;
end
AA_eq_g(1,:) = aux_r_g;
for ii=1:(TT-1)
    AA_eq_g(ii+1,:) = circshift(aux_r_g,ii);
end

%Construction of AA_eq for the STORAGES
aux_r_s = zeros(1,2*SS*TT);
aux_r_s(1,1)=1;
for ii=1:(SS-1)
    aux_r_s(1,(ii*(2*TT))+1)=1;
end
AA_eq_s(1,:) = aux_r_s;
for ii=1:(TT-1)
    AA_eq_s(ii+1,:) = circshift(aux_r_s,ii);
end

% Construction of AA_eq for the CONTROLLABLE_LOADs
aux_r_c = zeros(1,2*CC*TT);
aux_r_c(1,1)=1;
for ii=1:(CC-1)
    aux_r_c(1,(ii*(2*TT))+1)=1;
end
AA_eq_c(1,:) = aux_r_c;
for ii=1:(TT-1)
    AA_eq_c(ii+1,:) = circshift(aux_r_c,ii);
end

% Construction of AA_eq for the POWER TRADING COST (only last agent)
AA_eq_t = [eye(TT,TT),zeros(TT,TT)];

% Costruction of Final AA_eq and DD (that is the b for quadprog!)
AA_eq = [AA_eq_g,AA_eq_s,AA_eq_c,AA_eq_t]; %to be modified

time = [0:TT-1];
DD = (sin(time)+5)';

%Definition of the CENTRALIZED OPTIMIZATION
[zopt,fopt,exit_flag] = quadprog(2*QQ_centr,RR_centr,AA_tot,bb_tot,AA_eq,DD,[],[],[],options_qp);

if exit_flag ~= 1
  fprintf(2,'The problem %.4g problem occurred in the centralized solution\n',exit_flag);
  return;
end

fprintf('Centralized optimal cost is %.4g\n',fopt);

%Extraction of Popt from the augmented state variable Zopt 
popt = [];
popt(1:TT,1)=zopt(1:TT,1);
for ii=1:NN-1
    popt=[popt;zopt((ii*2*TT)+1:(ii*2*TT)+TT,1)];
end

%% Distributed solution

% Initializations
MAXITERS = 15*1e3;

LMD = rand(TT,NN,MAXITERS);
PP = zeros(TT,NN,MAXITERS);
dd = zeros(TT,NN,MAXITERS);

ZZ = zeros(2*TT,CC,MAXITERS); %variables for CONLs

KK = zeros(2*TT,MAXITERS); %variable for POWER TRADING COST

primal_cost = zeros(MAXITERS,1);
dual_cost   = zeros(MAXITERS,1);

consensus_err = zeros(MAXITERS,1);
tracking_err  = zeros(MAXITERS,1);
state_diff = zeros(MAXITERS,1);
feasibility = zeros(TT,MAXITERS);

LMD_init = rand(TT,NN);
LMD(:,:,1)= LMD_init;

%State initialization
for ii= 1:GG 
        p_init = quadprog(2*QQ_g(:,:,ii),RR_g(:,ii),AA_g(:,:,ii),bb_g(:,ii),[],[],[],[],[],options_qp);
        PP(:,ii,1) = p_init;
end
for ii= GG+1:SS+GG 
        p_init = quadprog(2*QQ_s(:,:,ii-GG),RR_s(:,ii-GG),AA_s(:,:,ii-GG),bb_s(:,ii-GG),[],[],[],[],[],options_qp);
        PP(:,ii,1) = p_init;
end
for ii= GG+SS+1:NN-1 
        ZZ(:,ii-GG-SS,1) = quadprog(2*QQ_c(:,:,ii-GG-SS),RR_c(:,ii-GG-SS),AA_c(:,:,ii-GG-SS),bb_c(:,ii-GG-SS),[],[],[],[],[],options_qp);
        p_init = ZZ(1:TT,ii-GG-SS,1);
        PP(:,ii,1) = p_init;
end
KK(:,1) = quadprog(2*QQ_t,RR_t,AA_t,bb_t,[],[],[],[],[],options_qp);
p_init = KK(1:TT,1);
PP(:,NN,1) = p_init;

%Gradient initialization
for ii = 1:NN-1
    dd(:,ii,1)=PP(:,ii,1); 
end
dd(:,NN,1)=PP(:,NN,1)-DD;

%Initialization of a constant stepsize
gamma = 6*1e-3;

for tt = 1:MAXITERS-1
  fprintf('iteration %d\n',tt);
  jj=0;
  %update of lamda
  for ii=1:NN
    N_ii = find(Adj(:,ii) == 1)';
    
     % LMD_i^t+1 =  sum_{j \in N_i U { i } } w_ij LMD_j^t + gamma dd_i^t
    LMD(:,ii,tt+1) = WW(ii,ii)*LMD(:,ii,tt); %the WW are what we called a_ij
    
    for jj = N_ii
      LMD(:,ii,tt+1) = LMD(:,ii,tt+1) + WW(ii,jj)*LMD(:,jj,tt); 
    end
    
    LMD(:,ii,tt+1) = LMD(:,ii,tt+1) + gamma*dd(:,ii,tt);
  end
  
  %update of PP, the state
  for ii = 1:GG 
         PP(:,ii,tt+1) = quadprog(2*QQ_g(:,:,ii),RR_g(:,ii)+LMD(:,ii,tt+1),AA_g(:,:,ii),bb_g(:,ii),[],[],[],[],[],options_qp);
  end
  for ii = GG+1:SS+GG 
         PP(:,ii,tt+1) = quadprog(2*QQ_s(:,:,ii-GG),RR_s(:,ii-GG)+LMD(:,ii,tt+1),AA_s(:,:,ii-GG),bb_s(:,ii-GG),[],[],[],[],[],options_qp);
  end
  for ii = GG+SS+1:NN-1 
         ZZ(:,ii-GG-SS,tt+1) = quadprog(2*QQ_c(:,:,ii-GG-SS),RR_c(:,ii-GG-SS)+[LMD(:,ii,tt+1);zeros(TT,1)],AA_c(:,:,ii-GG-SS),bb_c(:,ii-GG-SS),[],[],[],[],[],options_qp);
         PP(:,ii,tt+1) = ZZ(1:TT,ii-GG-SS,tt+1);
  end
  KK(:,tt+1) = quadprog(2*QQ_t,RR_t+[LMD(:,NN,tt+1);zeros(TT,1)],AA_t,bb_t,[],[],[],[],[],options_qp);
  PP(:,NN,tt+1) = KK(1:TT,tt+1);
  
  %update of dd, the estimate of the gradient
  for ii=1:NN
    N_ii = find(Adj(:,ii) == 1)';
    
     % dd_i^t+1 =  sum_{j \in N_i U { i } } w_ij dd_j^t + (p_i^t+1 -
     % p_i^t)
    dd(:,ii,tt+1) = WW(ii,ii)*dd(:,ii,tt); %the WW are what we called a_ij
    
    for jj = N_ii
      dd(:,ii,tt+1) = dd(:,ii,tt+1) + WW(ii,jj)*dd(:,jj,tt);
    end
    
    dd(:,ii,tt+1) = dd(:,ii,tt+1) + PP(:,ii,tt+1)- PP(:,ii,tt);
  end
  
  % Performance check
  LMD_avg = mean(LMD(:,:,tt),2);
  dd_avg = mean(dd(:,:,tt),2);
   for ii=1:GG  
    [ff_ii,~] = cost_functions(PP(:,ii,tt),QQ_g(:,:,ii),RR_g(:,ii));
    primal_cost(tt) = primal_cost(tt) + ff_ii;
    
    qq_ii = ff_ii + LMD(:,ii,tt)'*(PP(:,ii,tt));
    dual_cost(tt) = dual_cost(tt) + qq_ii;

    tracking_err(tt)  = tracking_err(tt)  + norm(dd(:,ii,tt) - dd_avg);
    consensus_err(tt) = consensus_err(tt) + norm(LMD(:,ii,tt) - LMD_avg);
    
    state_diff(tt) = state_diff(tt) + norm(PP(:,ii,tt)-popt((1+TT*(ii-1)):(TT*ii)));
    feasibility(:,tt) = feasibility(:,tt) + PP(:,ii,tt);
   end
  for ii=GG+1:SS+GG  
    [ff_ii,~] = cost_functions(PP(:,ii,tt),QQ_s(:,:,ii-GG),RR_s(:,ii-GG));
    primal_cost(tt) = primal_cost(tt) + ff_ii;
    
    qq_ii = ff_ii + LMD(:,ii,tt)'*(PP(:,ii,tt));
    dual_cost(tt) = dual_cost(tt) + qq_ii;

    tracking_err(tt)  = tracking_err(tt)  + norm(dd(:,ii,tt) - dd_avg);
    consensus_err(tt) = consensus_err(tt) + norm(LMD(:,ii,tt) - LMD_avg);
    
    state_diff(tt) = state_diff(tt) + norm(PP(:,ii,tt)-popt((1+TT*(ii-1)):(TT*ii)));
    feasibility(:,tt) = feasibility(:,tt) + PP(:,ii,tt);
  end
  for ii=GG+SS+1:NN-1  
    [ff_ii,~] = cost_functions(ZZ(:,ii-GG-SS,tt),QQ_c(:,:,ii-GG-SS),RR_c(:,ii-GG-SS));
    primal_cost(tt) = primal_cost(tt) + ff_ii;
    
    qq_ii = ff_ii + LMD(:,ii,tt)'*(PP(:,ii,tt));
    dual_cost(tt) = dual_cost(tt) + qq_ii;

    tracking_err(tt)  = tracking_err(tt)  + norm(dd(:,ii,tt) - dd_avg);
    consensus_err(tt) = consensus_err(tt) + norm(LMD(:,ii,tt) - LMD_avg);
    
    state_diff(tt) = state_diff(tt) + norm(PP(:,ii,tt)-popt((1+TT*(ii-1)):(TT*ii)));
    feasibility(:,tt) = feasibility(:,tt) + PP(:,ii,tt);
  end
   [ff_ii,~] = cost_functions(KK(:,tt),QQ_t,RR_t);
   primal_cost(tt) = primal_cost(tt) + ff_ii; 
   qq_ii = ff_ii + LMD(:,NN,tt)'*((PP(:,NN,tt))-DD);
   dual_cost(tt) = dual_cost(tt) + qq_ii;
   
   tracking_err(tt)  = tracking_err(tt)  + norm(dd(:,NN,tt) - dd_avg);
   consensus_err(tt) = consensus_err(tt) + norm(LMD(:,NN,tt) - LMD_avg);
   
   state_diff(tt) = state_diff(tt) + norm(PP(:,NN,tt)-popt((1+TT*(NN-1)):(TT*NN)));
   feasibility(:,tt) = feasibility(:,tt) + PP(:,NN,tt)-DD;
end

for ii=1:GG
    
    [ff_ii,~] = cost_functions(PP(:,ii,MAXITERS),QQ_g(:,:,ii),RR_g(:,ii));
 
    
    primal_cost(MAXITERS) = primal_cost(MAXITERS) + ff_ii;
    
    qq_ii = ff_ii + LMD(:,ii,MAXITERS)'*(PP(:,ii,MAXITERS));
    dual_cost(MAXITERS) = dual_cost(MAXITERS) + qq_ii;

    tracking_err(MAXITERS)  = tracking_err(MAXITERS)  + norm(dd(:,ii,MAXITERS) - dd_avg);
    consensus_err(MAXITERS) = consensus_err(MAXITERS) + norm(LMD(:,ii,MAXITERS) - LMD_avg);
    
    state_diff(MAXITERS) = state_diff(MAXITERS) + norm(PP(:,ii,MAXITERS)-popt((1+TT*(ii-1)):(TT*ii)));
    feasibility(:,MAXITERS) = feasibility(:,MAXITERS) + PP(:,ii,MAXITERS);
end
for ii=GG+1:SS+GG
    
    [ff_ii,~] = cost_functions(PP(:,ii,MAXITERS),QQ_s(:,:,ii-GG),RR_s(:,ii-GG));
 
    
    primal_cost(MAXITERS) = primal_cost(MAXITERS) + ff_ii;
    
    qq_ii = ff_ii + LMD(:,ii,MAXITERS)'*(PP(:,ii,MAXITERS));
    dual_cost(MAXITERS) = dual_cost(MAXITERS) + qq_ii;

    tracking_err(MAXITERS)  = tracking_err(MAXITERS)  + norm(dd(:,ii,MAXITERS) - dd_avg);
    consensus_err(MAXITERS) = consensus_err(MAXITERS) + norm(LMD(:,ii,MAXITERS) - LMD_avg);
    
    state_diff(MAXITERS) = state_diff(MAXITERS) + norm(PP(:,ii,MAXITERS)-popt((1+TT*(ii-1)):(TT*ii)));
    feasibility(:,MAXITERS) = feasibility(:,MAXITERS) + PP(:,ii,MAXITERS);
end
for ii=GG+SS+1:NN-1
    
    [ff_ii,~] = cost_functions(ZZ(:,ii-GG-SS,MAXITERS),QQ_c(:,:,ii-GG-SS),RR_c(:,ii-GG-SS));
 
    
    primal_cost(MAXITERS) = primal_cost(MAXITERS) + ff_ii;
    
    qq_ii = ff_ii + LMD(:,ii,MAXITERS)'*(PP(:,ii,MAXITERS));
    dual_cost(MAXITERS) = dual_cost(MAXITERS) + qq_ii;

    tracking_err(MAXITERS)  = tracking_err(MAXITERS)  + norm(dd(:,ii,MAXITERS) - dd_avg);
    consensus_err(MAXITERS) = consensus_err(MAXITERS) + norm(LMD(:,ii,MAXITERS) - LMD_avg);
    
    state_diff(MAXITERS) = state_diff(MAXITERS) + norm(PP(:,ii,MAXITERS)-popt((1+TT*(ii-1)):(TT*ii)));
    feasibility(:,MAXITERS) = feasibility(:,MAXITERS) + PP(:,ii,MAXITERS);
end
   [ff_ii,~] = cost_functions(KK(:,MAXITERS),QQ_t,RR_t); 
   primal_cost(MAXITERS) = primal_cost(MAXITERS) + ff_ii;
   qq_ii = ff_ii + LMD(:,NN,MAXITERS)'*((PP(:,NN,MAXITERS))-DD);
   dual_cost(MAXITERS) = dual_cost(MAXITERS) + qq_ii;
   
   tracking_err(MAXITERS)  = tracking_err(MAXITERS)  + norm(dd(:,NN,MAXITERS) - dd_avg);
   consensus_err(MAXITERS) = consensus_err(MAXITERS) + norm(LMD(:,NN,MAXITERS) - LMD_avg);
   
   state_diff(MAXITERS) = state_diff(MAXITERS) + norm(PP(:,NN,MAXITERS)-popt((1+TT*(NN-1)):(TT*NN)));
   
   feasibility(:,MAXITERS) = feasibility(:,MAXITERS) + PP(:,NN,MAXITERS)-DD;
   
   
   
    %%
figure
  semilogy(1:MAXITERS,abs(primal_cost(1:MAXITERS)-fopt), 'LineWidth',2);
  hold on, grid on, zoom on
  semilogy(1:MAXITERS,abs(dual_cost(1:MAXITERS)-fopt), 'LineWidth',2);
  xlabel('t')
  ylabel('cost error')
  legend('primal cost','dual cost')
  
  %%
figure
  semilogy(1:MAXITERS,state_diff(1:MAXITERS), 'LineWidth',2);
  hold on, grid on
  xlabel('t')
  ylabel('state error')
%%  
 figure
  semilogy(1:MAXITERS-1,consensus_err(1:MAXITERS-1), 'LineWidth',2);
  hold on, grid on, zoom on
  xlabel('t')
  ylabel('consensus error')
  
  %%
figure
  semilogy(1:MAXITERS,tracking_err(1:MAXITERS), 'LineWidth',2);
  hold on, grid on
  xlabel('t')
  ylabel('tracking error')
  
  
  %%
  feasibile_norm = [];
 for ii=1:MAXITERS
     feasibile_norm(ii) = norm(feasibility(:,ii));
 end
figure
  semilogy(1:MAXITERS,feasibile_norm(1:MAXITERS), 'LineWidth',2);
  hold on, grid on
  xlabel('t')
  ylabel('|| p_1 + ... + p_N - D ||')


