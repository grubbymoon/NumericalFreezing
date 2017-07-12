function J = calc_Jacobi(u_old, LHS, RHS, residual)

% set the increment
h = 1e-6;
% get the size of unknown vector
n = size(u_old,1);

for i = 1:n
    u_tmp = u_old;
    % evaluate the residual using 1st order forward Eular approximation
     u_tmp(i) = u_tmp(i) + h;
     residual_new = RHS - LHS*u_tmp;
     J(:,i) = (residual_new - residual)./h;
end
 
