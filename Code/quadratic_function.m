function [f_val,f_grad] = quadratic_function(x,Q,r)
  f_val  = x'*Q*x+r'*x;
  f_grad = 2*Q*x+r;
end 