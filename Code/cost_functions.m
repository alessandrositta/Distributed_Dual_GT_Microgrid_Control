function [f_val,f_grad] = cost_functions(p,Q,r)
%function only of the GENERATOR
  f_val  = p'*Q*p+r'*p;
  f_grad = (Q+Q')*p+r; %attention
end 