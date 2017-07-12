function [LHS,RHS] = LHS_RHS_Construction(T_ele, T_pre, theta, dt, coord)
global lambda_w;
global C_w;
global rho_w;
global lambda_i; 
global C_i;
global rho_i;
global lambda_s; 
global C_s;
global rho_s; 
global w;
global phi;
global L_I;

% Ice volume fraction with temperature
logistic = 1/(1 + exp( -w*T_ele));
phi_i = phi*(1 - logistic);
% determine the overall heat capacity in this element
C_p = rho_s*C_s*(1 - phi) + rho_i*C_i*phi_i + rho_w*C_w*(phi - phi_i) ;
% determine the overall thermal conductivity in this element
lambda = lambda_s*(1 - phi) + lambda_i*phi_i + lambda_w*(phi - phi_i) ;
% Ice volume fraction changing rate
% for logistic function y_dot = (1 - y)*y
sigmoid_dot =  - phi*w*logistic*(1 - logistic);


% local mass marix
M11 = (C_p - rho_i*sigmoid_dot*L_I);

% FEM
FEM_M = (1/dt) * M11 * shapeshape_tri(coord);

% ------------------------------------------------------------------------

% local transfer matrix 
D22 = lambda ;
% FEM
FEM_K = dshapedshape_tri(coord, [D22, 0; 0, D22]);


% assemble to the LHS and RHS
LHS = FEM_M+theta*FEM_K;
RHS = (FEM_M - ((1-theta) * FEM_K))* T_pre;  

 
end