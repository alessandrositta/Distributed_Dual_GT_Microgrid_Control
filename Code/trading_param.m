function [AA_tot,bb_tot,QQ_centr,RR_centr,AA_t,bb_t,QQ_t,RR_t] = trading_param(TT,epsilon,c1,c2,E,AA_tot_in,bb_tot_in,QQ_centr_in,RR_centr_in)
% POWER TRADING COST parameters initialization
% Given the function and gradient for the POWER TRADING COST:
% f_i(p_i) = p_i^T Q_i p_i + R_i^T p_i)
% grad_f_i(p_i) = (Q_i+Q_i^T) p_i + R_i)
% This function RETURNS the matrices (consider quadprog):
% AA_tot --> block diagonal matrix, describes the unequality constraint A 
% bb_tot --> column vector, describes the unequality constraint b 
% QQ_centr --> block diagonal matrix, describe the cost function matrix Q
% RR_centr --> matrix, describe the cost function matrix R
%
% Good values for the inputs are: 
% epsilon = 1
% c1 = 2
% c2 = 3
% E = 5

% Matrices QQ and RR of the cost function for the POWER TRADING COST
aux = zeros(2*TT,2*TT);
for ii=1:TT
    aux(ii,ii) = epsilon;
end

QQ_t = aux;
RR_t = [-c1*ones(TT,1);c2*ones(TT,1)];

% Constraints matrix creation 
AA_1_t = [eye(TT,TT),-eye(TT,TT)];
bb_1_t = zeros(TT,1); 

AA_2_t = -[zeros(TT,TT),eye(TT,TT)];
bb_2_t = zeros(TT,1);

AA_3_t = [eye(TT,TT),zeros(TT,TT)];
bb_3_t = E*ones(TT,1);

AA_4_t = -[eye(TT,TT),zeros(TT,TT)];
bb_4_t = E*ones(TT,1);

AA_tot = AA_tot_in;
bb_tot = bb_tot_in;
QQ_centr = QQ_centr_in;
RR_centr = RR_centr_in;

% Matrices containing the local constraints for POWER TRADING COST
AA_t = [AA_1_t;AA_2_t;AA_3_t;AA_4_t];
bb_t = [bb_1_t;bb_2_t;bb_3_t;bb_4_t];

% Matrices for the CENTRALIZED SOLUTION

  AA_tot = blkdiag(AA_tot,[AA_1_t;AA_2_t;AA_3_t;AA_4_t]);
  bb_tot = [bb_tot;bb_1_t;bb_2_t;bb_3_t;bb_4_t];
  QQ_centr = blkdiag(QQ_centr,QQ_t);
  RR_centr = [RR_centr;RR_t];


end 