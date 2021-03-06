% % % % % % % % % % % % % % % % % % % % % % % % % 
% freezing equation solver 1D (picards method for
% nonlinearlarity)  (Equilibrium Approach) sigmoid
%  
% Authors: Tianyuan Zheng
% Email  : tianyuan.zheng@ufz.de
% 
% % % % % % % % % % % % % % % % % % % % % % % % 

% The governing equation is: 
% C_p*dT/dt + L_I*rhophrate_I - lambda*d2T/dx2 = 0 
% rhophrate_I = alpha*(T - T_m) = rho_IR dphi_i/dt
% lambda = lambda_s*(1 - phi) + lambda_i*phi_i + lambda_w*(phi - phi_i)
% C_p = rho_s*C_s*(1 - phi) + rho_i*C_i*phi_i + rho_w*C_w*(phi - phi_i)


clear;
% parameter input ----------------------------
% define heat capacity of water
global C_w;
C_w = 4187 ;   % J*kg-1*K-1

% define heat capacity of ice 
global C_i;
C_i = 2108 ;

% define heat capacity of soil
global C_s;
C_s = 840 ;

% define thermal conductivity of water
global lambda_w;
lambda_w  = 0.58 ;  % W*m-1*K-1

% define thermal conductivity of ice
global lambda_i;
lambda_i = 2.14 ;

% define thermal conductivity of soil
global lambda_s;
lambda_s = 2.9 ;

% define density of  water 
global rho_w;
rho_w = 997 ; % kg*m3-1

% define density of ice 
global rho_i; 
rho_i = 997 ; % since no volume expansion considered, 

% define density of soil
global rho_s;
rho_s = 2600 ;

% define the porosity 
global phi;
phi = 0.05 ;

% define the latent heat of phase change
global L_I;
L_I = 334000 ; % J/kg

% define the melting point 
global T_m;
T_m = 0 ;% Celsius degree

% define the parameters for Sigmoid function
global w;
w = 5 ;

% Preprocessing -------------------------------
% ------------- reading the mesh ----- node and element ----------
format long;
nodes      = dlmread('2d_test_node.txt');
nodes      = nodes(:,:);

% reading the element connectivity
elements   = dlmread('2d_test_element.txt');
connectivity=dlmread('2d_test_element.txt');
% convert to matlab index start from one
elements   = elements + 1 ; 
% now the spacial dimension and element type
[nn,sdim]  = size(nodes) ;
[ne,dummy] = size(elements) ;

% ------------------ boundary conditions -----------------------
 % get the index of nodes that the x value equals 0 
index1_lhs = find(nodes(:,1) == 0);
index2_lhs=nn+index1_lhs;
index_lhs=[index1_lhs;index2_lhs];
% define this node to the boundary node
index_bcn = index1_lhs ;
% assign the x,y value of the boundary into 'bc_middle_nodes'
bc_lhs_nodes = nodes(index_bcn,:);
% assign the number of row and column of 'bc_middle_nodes' to 'nbcn' and 'dummy'
[nbcn,dummy] = size( bc_lhs_nodes );
bc_nodes_value = -6 * ones(nbcn,1);
% bc_ambient_nodes are the rest part
% index = find(bc_all_nodes ~= bc_middel_nodes);
% bc_ambient_nodes = bc_all_nodes(index,:);

% --------------------- Initial conditions -------------------------
% ------ T is for every node, other parameters are for elements ---
T_ini = ones(nn,1); 
T_ini = 4 * T_ini;
T_ini(index1_lhs,1) = -6 ;

% initialize previous and current solution matrix
U_pre_timestep  = T_ini;
U_pre_Picard  = T_ini;
U_cur_Picard = U_pre_Picard ;

% initialize parameter
phi_i = sparse(ne,1 );
C_p = sparse(ne,1 );
lambda = sparse(ne,1);
% at first no ice existed 
phi_i = sparse(ne,1);

% Processing ------------------------------------------------------
% time control
time_steps = [0:90:86400]';
[steps,dummy] = size(time_steps); 
steps = steps - 1; 
theta = 1.0;  % implicit
% theta = 0.5;  % Crank-Nicolson


% storage space of unknowns after each time step
U_record      = zeros(nn, steps+1);
U_record(:,1) = U_pre_Picard;
phi_i_record =  zeros(ne, steps+1);
phi_i_record(:,1) = phi_i;

% ------------- iteration control -----------------
tol = 1e-3 ; % set the tolerance of picard
counter = 0 ; % count the times of iteration
maxiter = 40 ;
max_lin_sol_iter = 1000; 
    
% time step loop 
for ti = 1 : steps
    
    str = ['Time step ', num2str(ti), ' at time ', num2str(time_steps(ti+1))];
    disp(str);
    dt = time_steps(ti+1) - time_steps(ti); 
    error = 100 ; % set the initial error
    counter = 0 ;
    
   while (error > tol) && (counter < maxiter)
    % cleaning of LHS and RHS
    LHS    = sparse(nn,nn);
    RHS    = sparse(nn,1 );
    % initial guess of the temperature of the element based on nodes
    U_pre_Picard = U_cur_Picard;
    
    % loop over all the elements, 
    for ie = 1 : ne % 'ne' is the number of the total elements 
        % get the coordinates of connecting nodes
        sctr  = elements(ie,:) ; % get the nodes of the local element 
        coord = nodes(sctr,:) ; % get the coordinates of the nodes 
        %Interpolation in a tri element
        [nrows,ncols] = size(coord);
        I = ones(nrows,1);
        A = 0.5 * det([I coord]);
        U_avg=(1/A)*[A/3,A/3,A/3]*U_pre_Picard(sctr);
        % assemble the local matrix 
        [LHS_loc,RHS_loc] = LHS_RHS_Construction(U_avg, U_pre_timestep(sctr), theta, dt, coord);
        % assemble the global matrix
        LHS(sctr,sctr) = LHS(sctr,sctr) + LHS_loc;  
        RHS(sctr) = RHS(sctr) + RHS_loc; 
    end


    
    % check out the sparsity of LHS
    % spy(LHS);
    % imposing boundary conditions---------------------
    % loop over all boundary condition nodes
    for ib = 1 : nbcn
        idx = index1_lhs(ib);
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
    
    % solve linear equation system. 
    U_cur_Picard = bicgstab(LHS,RHS,1e-6,max_lin_sol_iter) ; % solve the linear system
    % u_cur = LHS \ RHS;
    error_list = U_cur_Picard - U_pre_Picard ;
    error = norm(error_list) ;
    counter = counter + 1 ;
    disp(error);
   end
   
    % ---------------------- Convert to VTK format -----------------
    x=nodes(:,1);
    y=nodes(:,2);
    legen=['Temperature'];
    txt=['dataU',num2str(ti),'.vtk'];
    c_nodal=(full(U_cur_Picard(1:nn)))';
    matlab2vtk(x,y,connectivity,txt,c_nodal,1,nn,ne,legen);
    
   
    U_pre_timestep=U_cur_Picard;
    % save current result
    U_record(:,ti+1) = U_cur_Picard;  
    phi_i_record(:,ti+1) = phi_i;
    % output iteration times
    str2 = ['iteration times ', num2str(counter)];
    disp(str2);
    
    
    
    % hold on;
end

% Postprocessing-------------------------------
% 
% hold on;


% scatter(nodes(:,1),nodes(:,2),5,u_record(:,10));

% plot(nodes(:,1),u_record(:,1),nodes(:,1),u_record(:,151),nodes(:,1),u_record(:,301),nodes(:,1),u_record(:,451),nodes(:,1),u_record(:,601),nodes(:,1),u_record(:,1201),nodes(:,1),u_record(:,1801),nodes(:,1),u_record(:,2041),nodes(:,1),u_record(:,2101),nodes(:,1),u_record(:,2161),nodes(:,1),u_record(:,2221),nodes(:,1),u_record(:,2281),nodes(:,1),u_record(:,2341),nodes(:,1),u_record(:,2401),nodes(:,1),u_record(:,2461),nodes(:,1),u_record(:,2521),nodes(:,1),u_record(:,2581),nodes(:,1),u_record(:,2641),nodes(:,1),u_record(:,2701),nodes(:,1),u_record(:,2761),nodes(:,1),u_record(:,2821),nodes(:,1),u_record(:,2881),nodes(:,1),u_record(:,2941),nodes(:,1),u_record(:,3001))


% plot(nodes(:,1),u_record(:,1),nodes(:,1),u_record(:,1501),nodes(:,1),u_record(:,3001),nodes(:,1),u_record(:,4501),nodes(:,1),u_record(:,6001),nodes(:,1),u_record(:,6601),nodes(:,1),u_record(:,7201),nodes(:,1),u_record(:,7431),nodes(:,1),u_record(:,7491),nodes(:,1),u_record(:,7551),nodes(:,1),u_record(:,7611),nodes(:,1),u_record(:,7671),nodes(:,1),u_record(:,7731),nodes(:,1),u_record(:,7791),nodes(:,1),u_record(:,7851),nodes(:,1),u_record(:,7911),nodes(:,1),u_record(:,7971),nodes(:,1),u_record(:,8031),nodes(:,1),u_record(:,8091),nodes(:,1),u_record(:,8151),nodes(:,1),u_record(:,8211),nodes(:,1),u_record(:,8271),nodes(:,1),u_record(:,8331),nodes(:,1),u_record(:,8391))





