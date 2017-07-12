 % % % % % % % % % % % % % % % % % % % % % % % % 
% TH coupled freezing equation solver (picards method for
% nonlinearlarity) (Equilibrium Approach) sigmoid (Assume K and Ss are constant) 
%
% 
% Authors: Tianyuan Zheng
% Email  : tianyuan.zheng@ufz.de
% 
% % % % % % % % % % % % % % % % % % % % % % % % 

% The governing equation is: 
% For the thermal part:
% C_p*dT/dt + L_I*rhophrate_I - lambda*d2T/dx2 + rho_f*C_f*(-K*dh/dx)= 0 
% rhophrate_I = alpha*(T - T_m) = rho_IR dphi_i/dt
% lambda = lambda_s*(1 - phi) + lambda_i*phi_i + lambda_w*(phi - phi_i)
% C_p = rho_s*C_s*(1 - phi) + rho_i*C_i*phi_i + rho_w*C_w*(phi - phi_i)
% For the hydraulic part:
% Ss*dh/dt - K*d2h/dh2 = 0
% Ss = Constant 
% K = 


% parameter input ----------------------------
% define heat capacity of water
global C_w;
C_w = 4187 ;   % J*kg-1*K-1
% define thermal conductivity of water
global lambda_w;
lambda_w  = 0.58 ;  % W*m-1*K-1
% define density of  water 
global rho_w;
rho_w = 997 ; % kg*m3-1
% define heat capacity of ice 
global C_i;
C_i = 2108 ;
% define thermal conductivity of ice
global lambda_i;
lambda_i = 2.14 ;
% define density of ice
global rho_i;
rho_i = 997 ; 
% define real density of ice (since we neglect volume change and
% mechanics of ice) we set the density of ice and water same)
global rho_IR;
rho_IR = 997 ;
% define heat capacity of soil
global C_s;
C_s = 840 ;
% define thermal conductivity of soil
global lambda_s; 
lambda_s = 2.9 ;
% define density of soil
global rho_s;
rho_s = 2600 ;
% define the latent heat of phase change
global L_I;
L_I = 334000 ; % J/kg
% define the melting point 
global T_m;
T_m = 0 ;% Celsius degree
% define the porosity
global phi;
phi = 0.5 ;
% define the parameters for Sigmoid function
global w;
w = 10 ;
% define the hydraulic conductivity
global K0; % function of porosity
K0 = 10e-6 ; % m/s
% define the Storage coefficient
global Ss; 
Ss = 3*10e-7 ; % function of porosity

% boundary temperature
%T_B_amb = -10;     % Celsius degree
%T_A_in  = 10;     % Celsius degree
% end of parameter input ----------------------

% Preprocessing -------------------------------
% reading the nodes
format long;
%nodes      = dlmread('nodes 1000.txt');
nodes      = dlmread('nodes_coordinates_line.txt');
nodes      = nodes(:,:);
% reading the element connectivity
%elements   = dlmread('elements 1000.txt');
elements   = dlmread('elem_connectivity_line.txt');
connectivity   = dlmread('elem_connectivity_line.txt');
% convert to matlab index start from one
elements   = elements + 1; 
% now the spacial dimension and element type
[nn,sdim]  = size(nodes);
[ne,dummy] = size(elements);


% -------------------------------------
%|                                     |
%|                                     |
%|                                     |
%|                                     |
% -------------------------------------
%find boundary node (1D)
%Left hand side
index1_lhs = find(nodes(:,1) == 0); % get the index of nodes that the x value equals 0
index2_lhs=index1_lhs+nn;
index_lhs=[index1_lhs;index2_lhs];
bc_lhs_nodes=nodes(index1_lhs,:);
[nbcn_lhs,dummy] = size( bc_lhs_nodes );
%index2 = find(nodes(:,2) > 0.075); % get the index of nodes that the y value bigger than 0.075
%index3 = find(nodes(:,2) < 0.125); % get the index of nodes that the y value smaller than 0.125
% index of the three regions above 
% index_bcn  = intersect(intersect(index1, index2), index3); 
bc_nodes_value_T = -6 * ones(nbcn_lhs,1);
bc_nodes_value_H = 11 * ones(nbcn_lhs,1);
bc_nodes_value = Vector_Combine(bc_nodes_value_T,bc_nodes_value_H);
% bc_ambient_nodes are the rest part
% index = find(bc_all_nodes ~= bc_middel_nodes);
% bc_ambient_nodes = bc_all_nodes(index,:);

% initial conditions (note! T and h is for every node, other parameters are for elements)
T_ini = ones(nn,1); 
T_ini = 4 * T_ini; 
H_ini = ones(nn,1);
H_ini = 10 * H_ini ;
U_ini = Vector_Combine(T_ini,H_ini) ;
phi_i_initial = 0*ones(ne,1); 

% Processing ----------------------------------
% numerical control
max_lin_sol_iter = 1000; 
time_steps = [0:900:86400]';
[steps,dummy] = size(time_steps); 
steps = steps - 1; 
theta = 1.0;  % implicit
% theta = 0.5;  % Crank-Nicolson

% initialize the previous time step value and previous iteration step value
% Primary unknown ------T means temperature
% Primary unknown ------H means water head
T_pre_time  = T_ini;
T_pre_Picard  = T_ini;
H_pre_time  = H_ini;
H_pre_Picard  = H_ini; 
U_pre_time= U_ini;
U_pre_Picard= U_ini;
U_cur_Picard= U_pre_Picard;

% impose the Boundary condition 
for ib=1:2*nbcn_lhs
    idx=index_lhs(ib);
    U_pre_Picard(idx)=bc_node_value_U(ib);
    U_pre_time(idx)=bc_node_value_U(ib);
end

% initialize parameters
T_ele = sparse(ne,1 );
phi_i = sparse(ne,1 );
C_p = sparse(ne,1 );
lambda = sparse(ne,1);
K = sparse(ne,1);

tol = 1e-3 ; % set the tolerance of picard
counter = 0 ; % count the times of iteration
maxiter = 40 ;

% storage space of unknowns after each time step
T_record      = sparse(nn, steps+1);

H_record      = sparse(nn, steps+1);

% THE COUPLED UNKNOWN
U_RECORD=[T_record',H_record']';
Res=ones(2*nn,1);
    
% time step loop 
for ti = 1 : steps
str = ['Time step ', num2str(ti), ' at time ', num2str(time_steps(ti+1))];
disp(str);
dt = time_steps(ti+1) - time_steps(ti); 
Res = ones(2*nn,1); % set the initial residual
counter = 0 ;
    
   while (norm(Res) > tol) && (counter < maxiter)
   str=[num2str(norm(Res))];
   disp(str);
   U_pre_Picard=U_cur_Picard;
   T_pre_Picard=U_pre_Picard(1:nn);
   H_pre_Picard=U_pre_Picard(nn+1:2*nn);
   % cleaning of LHS and RHS
    LHS=sparse(2*nn,2*nn);
    RHS=sparse(2*nn,1);
    LHS11=sparse(nn,nn);
    LHS21= sparse(nn,nn);
    LHS12=sparse(nn,nn);
    LHS22=sparse(nn,nn);
    RHS_U=sparse(nn,1 );
    RHS_D=sparse(nn,1 );

    % loop over all the elements, 
    for ie = 1 : ne % 'ne' is the number of the total elements 
        % get the coordinates of connecting nodes
        sctr  = elements(ie,:) ; % get the nodes of the local element 
        coord = nodes(sctr,:) ; % get the coordinates of the nodes 
        L = coord(2) - coord(1) ; % get the length of the current element
        % initial guess of the temperature of the element based on nodes
        T_ele(ie) = mean(T_pre_Picard(sctr));
        H_ele(ie) = mean(H_pre_Picard(sctr));
        [LHS_loc,RHS_loc] = LHS_RHS_Construction(T_ele(ie), T_pre_time(sctr), H_pre_time(sctr), dt, theta, L);
        % Assemble the global LHS and RHS
        LHS11(sctr,sctr)= LHS11(sctr,sctr) + LHS_loc(1:2,1:2);
        LHS21(sctr,sctr)= LHS21(sctr,sctr) + LHS_loc(3:4,1:2);
        LHS12(sctr,sctr) = LHS12(sctr,sctr) + LHS_loc(1:2,3:4);
        LHS22(sctr,sctr) =LHS22(sctr,sctr) + LHS_loc(3:4,3:4);
        RHS_U(sctr) = RHS_U(sctr) + RHS_loc(1:2);
        RHS_D(sctr) = RHS_D(sctr) + RHS_loc(3:4);
    end
    %assembly each entry of the matrix
    LHS_COMP={LHS11 LHS12;LHS21 LHS22};
    LHS=cell2mat(LHS_COMP);
    RHS=Vector_Combine(RHS_U,RHS_D);
    RHS_dense=full(RHS);    

    % imposing boundary conditions---------------------
    % loop over all boundary condition nodes
    for ib = 1 : 2*nbcn_lhs
        idx = index_lhs(ib);
        % b(i) -= A(i,ii)*u_bar
        RHS = RHS - LHS(:,idx) * bc_nodes_value(ib);
        % A(ii,ii) -> xii
        xii = LHS(idx,idx);
        % A(ii,j ) =  0
        LHS(idx,:) = 0.0;
        % A(i ,ii) = 0
        LHS(:,idx) = 0.0;
        % b(ii) = xii * u_bar
        RHS(idx) = xii * bc_nodes_value(ib);
        % A(ii,ii) <- xii
        LHS(idx,idx) = xii; 
    end
    % end of imposing boundary conditions--------------
%     
%     [i,j,val] = find(LHS);
%     data_dump = [i,j,val];
%     dlmwrite('LHS_right.txt',data_dump);
%     [i,j,val] = find(RHS);
%     data_dump = [i,j,val];
%     dlmwrite('RHS_right.txt',data_dump);
    
    Res = RHS - LHS*U_pre_Picard ;
    % solve linear equation system. 
    U_cur_Picard = bicgstab(LHS,RHS,1e-6,max_lin_sol_iter) ; % solve the linear system
    % u_cur = LHS \ RHS;
    Res = U_cur_Picard - U_pre_Picard ;
    counter = counter + 1 ;
   end
    % show the iteration times in the last time step (have problems with 1D mesh)
    str2 = ['iteration times ', num2str(counter)];
    disp(str2);
    % output vtk file
    x=nodes(:,1);
    legen=['X'];
    txt=['dataX',num2str(ti),'.vtk'];
    c_nodal=(full(U_cur_Picard(1:nn)))';
    c_nodal2=(full(U_cur_Picard(nn+1:2*nn)))';
    matlab2vtk(x,connectivity,'data.vtk',c_nodal,1,nn,ne,legen);
    matlab2vtk(x,connectivity,txt,c_nodal2,1,nn,ne,legen);
    
    
    % hold on;
end

% Postprocessing-------------------------------
% 
% hold on;





