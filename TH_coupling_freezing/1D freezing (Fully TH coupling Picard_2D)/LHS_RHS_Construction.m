function [LHS,RHS] = LHS_RHS_Construction(T_ele, T_pre, H_pre, dt, theta, coord)
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
global K0 ;
global Ss ;

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
% Permeability change with porosity
K = ((phi-phi_i)/phi)^3*(1-phi)/(1-phi+phi_i)*K0 ;
% Get the gradient of H
%h_gradient = zeros(1,3) ;
%for i=1:3
%h_gradient(i) = dshapef_tri(coord, [H_pre(i),0;0,H_pre(i)]) ;
%end
%H_g =(1/A)*[A/3,A/3,A/3]*h_gradient;
H_g = 0.0001 ;

% local mass marix
M11 = (C_p - rho_i*sigmoid_dot*L_I);
M12 = 0;
M21 = 0;
M22 = Ss ;
% FEM
Fem_M11 = (1/dt)*M11*shapeshape_tri(coord);
Fem_M12 = (1/dt)*M12*shapeshape_tri(coord);
Fem_M21 = (1/dt)*M21*shapeshape_tri(coord);
Fem_M22 = (1/dt)*M22*shapeshape_tri(coord);
% -----------------------Eq 1

% local transfer matrix 
D11 = K ; 
D12 = 0 ;
D21 = 0 ;
D221 = -rho_w*C_w*K*H_g;
D222 = lambda ;
% FEM
Fem_K11 = dshapedshape_tri(coord, [D11, 0;0, D11]) ;
Fem_K12 = dshapedshape_tri(coord, [D12, 0;0, D12]) ;
Fem_K21 = dshapedshape_tri(coord, [D21, 0;0, D21]) ;
Fem_K22 = shapedshape_tri(coord, [D221, 0;0, D221]) + dshapedshape_tri(coord, [D222, 0;0, D222]) ;

%------------------------Eq 2
% assemble to the LHS
LHS11=Fem_M11+theta*Fem_K11;
LHS12=Fem_M12+theta*Fem_K12;
LHS21=Fem_M21+theta*Fem_K21;
LHS22=Fem_M22+theta*Fem_K22;

RHS1_U=(Fem_M12-(1-theta)*Fem_K12)*H_pre+(Fem_M11-(1-theta)*Fem_K11)*T_pre;
RHS1_D=(Fem_M22-(1-theta)*Fem_K22)*H_pre+(Fem_M21-(1-theta)*Fem_K21)*T_pre;
LHS_COMP={LHS11 LHS12;LHS21 LHS22};
LHS=cell2mat(LHS_COMP);
RHS=Vector_Combine(RHS1_U,RHS1_D);
 
end