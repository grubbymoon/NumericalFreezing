function J = calc_Jacobi(u_cur, LHS_j, RHS_j, LHS, RHS, residual, h, u_norm)


% get the size of unknown vector
n = size(u_cur,1);


switch 2
    case 1    
for i = 1:n
% evaluate the residual using 1st order forward Eular approximation
     u_tmp = u_cur;
if   u_norm > 1e-10;
     u_tmp(i) = u_tmp(i) + h*u_norm;
     residual_new = RHS_j - LHS_j*u_tmp;
     J(:,i) = (residual_new - residual)./(h*u_norm);
else 
     u_tmp(i) = u_tmp(i) + h;
     residual_new = RHS_j - LHS_j*u_tmp;
     J(:,i) = (residual_new - residual)./h;
end
end
    case 2
for i = 1:n
     u_tmp = u_cur;
% evaluate the residual using 1st order forward Eular approximation
if   u_norm > 1e-10;
     u_tmp(i) = u_tmp(i) + h*u_norm;
     residual_new = RHS - LHS*u_tmp;
     J(:,i) = (residual_new - residual)./(h*u_norm);
else 
     u_tmp(i) = u_tmp(i) + h;
     residual_new = RHS - LHS*u_tmp;
     J(:,i) = (residual_new - residual)./h;
end
end
end