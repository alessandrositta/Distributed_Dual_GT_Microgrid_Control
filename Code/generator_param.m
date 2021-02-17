function [AA_tot,bb_tot,QQ_centr,RR_centr,AA_g,bb_g,QQ_g,RR_g] = generator_param(TT,GG,alpha1,alpha2,ub_1_g,lb_2_g,ub_3_g,lb_4_g,AA_tot_in,bb_tot_in,QQ_centr_in,RR_centr_in)
% GENERATORS parameters initialization
% Given the function and gradient for the Generators:
% f_i(p_i) = p_i^T Q_i p_i + R_i^T p_i)
% grad_f_i(p_i) = (Q_i+Q_i^T) p_i + R_i)
% This function RETURNS the matrices (consider quadprog):
% AA_tot --> block diagonal matrix, describes the unequality constraint A 
% bb_tot --> column vector, describes the unequality constraint b 
% QQ_centr --> block diagonal matrix, describe the cost function matrix Q
% RR_centr --> matrix, describe the cost function matrix R
%
% Good values for the inputs are: 
% alpha1 = 1
% alpha1 = 2
% ub_1_g = 5
% lb_2_g = -1
% ub_3_g = 5
% lb_4_g = -1

% Matrices RR and QQ for the cost function of the GENERATORS
RR_g = alpha1*ones(TT,GG);
QQ_g = repmat(alpha2*eye(TT,TT), 1,1,GG);

% Constraints matrix creation 
AA_1_g = repmat(eye(TT,TT), 1,1,GG);
bb_1_g = ub_1_g*ones(TT,GG);

AA_2_g = repmat(-1*eye(TT,TT), 1,1,GG);
bb_2_g = -1*lb_2_g*ones(TT,GG);

Aux = -1*eye(TT-1,TT) + (triu(ones(TT-1,TT),1)-triu(ones(TT-1,TT),2));
AA_3_g = repmat(Aux, 1,1,GG);
bb_3_g = ub_3_g*ones(TT-1,GG);

AA_4_g = repmat(-Aux, 1,1,GG);
bb_4_g = -1*lb_4_g*ones(TT-1,GG);

AA_tot = AA_tot_in;
bb_tot = bb_tot_in;
QQ_centr = QQ_centr_in;
RR_centr = RR_centr_in;

% Matrices containing the local constraints for GENERATORS
AA_g = [AA_1_g;AA_2_g;AA_3_g;AA_4_g];
bb_g = [bb_1_g;bb_2_g;bb_3_g;bb_4_g];

% Matrices for the CENTRALIZED SOLUTION
for i = 1:GG
  AA_tot = blkdiag(AA_tot,[blkdiag(AA_1_g(:,:,i),zeros(TT,TT));blkdiag(AA_2_g(:,:,i),zeros(TT,TT));blkdiag(AA_3_g(:,:,i),zeros(TT-1,TT));blkdiag(AA_4_g(:,:,i),zeros(TT-1,TT))]);
  bb_tot = [bb_tot;bb_1_g(:,i);zeros(TT,1);bb_2_g(:,i);zeros(TT,1);bb_3_g(:,i);zeros(TT-1,1);bb_4_g(:,i);zeros(TT-1,1)];
  QQ_centr = blkdiag(QQ_centr,QQ_g(:,:,i),zeros(TT,TT)); %we added 0s to consider also the Y of the CONL
  RR_centr = [RR_centr;RR_g(:,i);zeros(TT,1)]; %we added 0s to consider also the Y of the CONL
end

end 
