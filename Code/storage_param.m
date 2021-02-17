function [AA_tot,bb_tot,QQ_centr,RR_centr,AA_s,bb_s,QQ_s,RR_s] = storage_param(TT,SS,epsilon,ub_1_s,lb_2_s,ub_3_s,lb_4_s,AA_tot_in,bb_tot_in,QQ_centr_in,RR_centr_in)
% STORAGE parameters initialization
% Given the function and gradient for the Storage:
% f_i(p_i) = p_i^T Q_i p_i 
% grad_f_i(p_i) = (Q_i+Q_i^T) p_i
% This function RETURNS the matrices (consider quadprog):
% AA_tot --> block diagonal matrix, describes the unequality constraint A 
% bb_tot --> column vector, describes the unequality constraint b 
% QQ_centr --> block diagonal matrix, describe the cost function matrix Q
% RR_centr --> matrix, describe the cost function matrix R
%
% Good values for the inputs are: 
% epsilon = 1
% ub_1_s = 3
% lb_2_s = -0.5
% ub_3_s = 4
% lb_4_s = 0

% Matrices QQ and RR of the cost function for the STORAGS
QQ_s = repmat(epsilon*eye(TT,TT), 1,1,SS);
RR_s = zeros(TT,SS);

% Constraints matrix creation 
AA_1_s = repmat(eye(TT,TT), 1,1,SS);
bb_1_s = ub_1_s*ones(TT,SS);

AA_2_s = repmat(-1*eye(TT,TT), 1,1,SS);
bb_2_s = -lb_2_s*ones(TT,SS);

Aux = -1*eye(TT-1,TT) + (triu(ones(TT-1,TT),1)-triu(ones(TT-1,TT),2));
Aux2 = pinv(Aux);
Aux3 = [eye(TT-1,TT-1),zeros((TT-1),1)];
AA_3_s = repmat(Aux2*Aux3, 1,1,SS); 
bb_3_s = ub_3_s*ones(TT,SS); 

AA_4_s = repmat(-Aux2*Aux3, 1,1,SS); 
bb_4_s = lb_4_s*ones(TT,SS); 

AA_tot = AA_tot_in;
bb_tot = bb_tot_in;
QQ_centr = QQ_centr_in;
RR_centr = RR_centr_in;

% Matrices containing the local constraints for STORAGES
AA_s = [AA_1_s;AA_2_s;AA_3_s;AA_4_s];
bb_s = [bb_1_s;bb_2_s;bb_3_s;bb_4_s];

% Matrices for the CENTRALIZED SOLUTION
for i = 1:SS %we added 0s to consider also the Y of the CONL
  AA_tot = blkdiag(AA_tot,[blkdiag(AA_1_s(:,:,i),zeros(TT,TT));blkdiag(AA_2_s(:,:,i),zeros(TT,TT));blkdiag(AA_3_s(:,:,i),zeros(TT,TT));blkdiag(AA_4_s(:,:,i),zeros(TT,TT))]);
  bb_tot = [bb_tot;bb_1_s(:,i);zeros(TT,1);bb_2_s(:,i);zeros(TT,1);bb_3_s(:,i);zeros(TT,1);bb_4_s(:,i);zeros(TT,1)];
  QQ_centr = blkdiag(QQ_centr,QQ_s(:,:,i),zeros(TT,TT));
  RR_centr = [RR_centr;RR_s(:,i);zeros(TT,1)];
end

end 
