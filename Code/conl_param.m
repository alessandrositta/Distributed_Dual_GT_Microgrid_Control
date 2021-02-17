function [AA_tot,bb_tot,QQ_centr,RR_centr,AA_c,bb_c,QQ_c,RR_c] = conl_param(TT,CC,epsilon,beta,P_des,P_bound,AA_tot_in,bb_tot_in,QQ_centr_in,RR_centr_in)
% CONTROLLABE LOAD parameters initialization
% Given the function and gradient for the CONLs:
% f_i(p_i) = p_i^T Q_i p_i + R_i^T p_i)
% grad_f_i(p_i) = (Q_i+Q_i^T) p_i + R_i)
% This function RETURNS the matrices (consider quadprog):
% AA_tot --> block diagonal matrix, describes the unequality constraint A 
% bb_tot --> column vector, describes the unequality constraint b 
% QQ_centr --> block diagonal matrix, describe the cost function matrix Q
% RR_centr --> matrix, describe the cost function matrix R
%
% Good values for the inputs are: 
% epsilon = 2
% beta = 3
% P_bound = 5


% Matrices QQ and RR of the cost function for the CONTROLLABLE_LOADs
aux = zeros(2*TT,2*TT);
for ii=1:TT
    aux(ii,ii) = epsilon;
end

QQ_c = repmat(aux, 1,1,CC);
RR_c = repmat([zeros(TT,1);beta*ones(TT,1)],1,CC);

% Constraints matrix creation 
AA_1_c = repmat(-[eye(TT,TT),eye(TT,TT)], 1,1,CC);
bb_1_c = P_des; %P_des must be of dimension (TT)*CC and it must be given in input

AA_2_c = repmat(-[zeros(TT,TT),eye(TT,TT)], 1,1,CC);
bb_2_c = zeros(TT,CC);

AA_3_c = repmat([eye(TT,TT),zeros(TT,TT)], 1,1,CC);
bb_3_c = P_bound*ones(TT,CC);

AA_4_c = repmat(-[eye(TT,TT),zeros(TT,TT)], 1,1,CC);
bb_4_c = P_bound*ones(TT,CC);

AA_tot = AA_tot_in;
bb_tot = bb_tot_in;
QQ_centr = QQ_centr_in;
RR_centr = RR_centr_in;

% Matrices containing the local constraints for CONTROLLABLE_LOADs
AA_c = [AA_1_c;AA_2_c;AA_3_c;AA_4_c];
bb_c = [bb_1_c;bb_2_c;bb_3_c;bb_4_c];

% Matrices for the CENTRALIZED SOLUTION
for i = 1:CC
  AA_tot = blkdiag(AA_tot,[AA_1_c(:,:,i);AA_2_c(:,:,i);AA_3_c(:,:,i);AA_4_c(:,:,i)]);
  bb_tot = [bb_tot;bb_1_c(:,i);bb_2_c(:,i);bb_3_c(:,i);bb_4_c(:,i)];
  QQ_centr = blkdiag(QQ_centr,QQ_c(:,:,i));
  RR_centr = [RR_centr;RR_c(:,i)];
end

end 