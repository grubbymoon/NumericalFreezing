 % % % % % % % % % % % % % % % % % % % % % % % % 
% freezing equation solver without mechanics (picards method for )
% 
% Authors: Tianyuan Zheng
% Email  : tianyuan.zheng@ufz.de
% 
% % % % % % % % % % % % % % % % % % % % % % % % 

% The governing equation is: 
% C_p*dT/dt + rho_w*L_I*rhophrate_I - lambda*(d2T/dy2 + d2T/dx2) = 0 
% rhophrate_i = alpha*(T - T_m) = rho_IR dphi_i/dt
% lambda = lambda_s*(1 - phi) + lambda_i*phi_i + lambda_w*(phi - phi_i)
% C_p = rho_s*C_s*(1 - phi) + rho_i*C_i*phi_i + rho_w*C_w*(phi - phi_i)



% parameter input ----------------------------
% define heat capacity of water
C_w = 4187 ;   % J*kg-1*K-1
% define thermal conductivity of water
lambda_w  = 0.58 ;  % W*m-1*K-1
% define density of water 
rho_w = 1000 ; % kg*m3-1
% define heat capacity of ice 
C_i = 2108 ;
% define thermal conductivity of ice
lambda_i = 2.14 ;
% define density of ice
rho_i = 900 ; 
% define real density of ice
rho_IR = 900 ;
% define heat capacity of soil
C_s = 840 ;
% define thermal conductivity of soil
lambda_s = 2.9 ;
% define density of soil
rho_s = 2600 ;
% define the porosity 
% n = 0.34 ; % without unit
% define the latent heat of phase change
L_I = 334000 ; % J/kg
% define the melting point 
T_m = 0 ;% Celsius degree
% define the porosity 
phi = 0.3 ;
% define the alpha 
% alpha = -1e-3;
alpha = -1e-3;
% define the overall heat capacity in one element 
% define the overall thermal conductivity in one element

% boundary temperature
%T_B_amb = -10;     % Celsius degree
%T_A_in  = 10;     % Celsius degree
% end of parameter input ----------------------

% Preprocessing -------------------------------
% reading the nodes
nodes      = dlmread('nodes_coordinates.txt');
% reading the element connectivity
elements   = dlmread('elem_connectivity.txt');
% convert to matlab index start from one
elements   = elements + 1; 
% now the spacial dimension and element type
[nn,sdim]  = size(nodes);
[ne,dummy] = size(elements);
% plot the mesh to check

% boundary conditions
index1 = find(nodes(:,1) == 0); % get the index of nodes that the x value equals 0
index2 = find(nodes(:,2) > 0.075); % get the index of nodes that the y value bigger than 0.075
index3 = find(nodes(:,2) < 0.125); % get the index of nodes that the y value smaller than 0.125
% index of the three regions above 
index_bcn  = intersect(intersect(index1, index2), index3); 
% assign the x,y value of the boundary into 'bc_middle_nodes'
bc_middel_nodes = nodes(index_bcn,:);
% assign the number of row and column of 'bc_middle_nodes' to 'nbcn' and 'dummy'
[nbcn,dummy] = size( bc_middel_nodes );
bc_nodes_value = -20 * ones(nbcn,1);
% bc_ambient_nodes are the rest part
% index = find(bc_all_nodes ~= bc_middel_nodes);
% bc_ambient_nodes = bc_all_nodes(index,:);

% initial conditions (note! T is for every node, other parameters are for elements)
T_ini = ones(nn,1); 
T_ini = 10.0 * T_ini; 
phi_i_initial = 0*ones(ne,1); 
C_p_initial = (C_w*phi + C_s*(1-phi))*ones(ne,1);
lambda_initial = (lambda_w*phi + lambda_s*(1-phi))*ones(ne,1);

%phi_w_inital = 0.3*ones(ne,1);


% Processing ----------------------------------
% numerical control
max_lin_sol_iter = 1000; 
% time control
% time_steps = [0:10:400, 500:100:10000, 11000:1000:86400]';
% time_steps = [0:10:400, 500:100:10000]';
time_steps = [0:50:10000]';
% time_steps = [0:864:8640]'; % time interval is 864s for every step
% time_steps = [0:10:100]'
% time_Steps = [0:864:86400]';
[steps,dummy] = size(time_steps); 
steps = steps - 1; 
theta = 1.0;  % implicit
% theta = 0.5;  % Crank-Nicolson

% initialize previous and current solution matrix
u_pre  = sparse(nn,1 ); 
u_cur  = sparse(nn,1 ); 
u_pre  = T_ini;
% initial guess 
u_cur  = T_ini;

% initialize parameters
T_ele = sparse(ne,1 );
phi_i = sparse(ne,1 );
C_p = sparse(ne,1 );
lambda = sparse(ne,1);
%phi_w = sparse(ne,1 );
phi_i = phi_i_initial;
phi_i_old = phi_i_initial;

C_p = C_p_initial;
lambda = lambda_initial;
%phi_w = phi_w_inital;
tol = 1e-6 ; % set the tolerance of picard
counter = 0 ; % count the times of iteration
maxiter = 100 ;

% storage space of unknowns after each time step
u_record      = zeros(nn, steps+1);
u_record(:,1) = u_pre;
phi_i_record =  zeros(ne, steps+1);
phi_i_record(:,1) = phi_i_initial;
C_p_record =    zeros(ne, steps+1);
C_p_record(:,1) = C_p_initial;
lambda_record = zeros(ne, steps+1);
lambda_record(:,1) = lambda_initial;
T_ele_record = zeros(ne, steps+1) ;
for ts = 1 : ne % change the initial node temperature to initial element temperature
    noele = elements(ts,:);
    T_ele_record(ts,1) = mean(u_cur(noele));
end
    
    

% time step loop 
for ti = 1 : steps
    % cleaning of LHS and RHS
    LHS    = sparse(nn,nn);
    RHS    = sparse(nn,1 );
    phi_i = sparse(ne,1);
    T_ele = sparse(ne,1);
    C_p = sparse(ne,1 );
    lambda = sparse(ne,1);
    % cleaning of phi
    str = ['Time step ', num2str(ti), ' at time ', num2str(time_steps(ti+1))];
    disp(str);
    dt = time_steps(ti+1) - time_steps(ti); 
    error = 100 ; % set the initial error
    counter = 0 ;
    
   while (error > tol) && (counter < maxiter)
    % loop over all the elements, 
    for ie = 1 : ne % 'ne' is the number of the total elements 
        % get the coordinates of connecting nodes
        sctr  = elements(ie,:) ; % get the nodes of the local element 
        coord = nodes(sctr,:) ; % get the coordinates of the nodes 
        % initial guess of the temperature of the element based on nodes
        u_comparion = u_cur;
        T_ele(ie) = mean(u_cur(sctr));
        % evaluate the freezing value
        % rhophrate_i = freezing(T, phi_i(ie));
        if phi_i_old(ie) <= 0 && T_ele(ie) > T_m ;
            rhophrate_i = 0; % when there is no ice, no more melting
        elseif phi - phi_i_old(ie) <= 0 && T_ele(ie) < T_m ;
            rhophrate_i = 0; % when there is no water, no more freezing
        else 
            rhophrate_i = alpha*(T_ele(ie) - T_m) ;
        end
        % evaluate the change of porosity
        phi_i(ie) = phi_i_old (ie) + dt*rhophrate_i/rho_IR;
        % make a limitation of the ice volume fraction
        phi_i(ie) = min(phi_i(ie),phi);
        phi_i(ie) = max(phi_i(ie),0);
        % determine the overall heat capacity in this element
        C_p(ie) = rho_s*C_s*(1 - phi) + rho_i*C_i*phi_i(ie) + rho_w*C_w*(phi - phi_i(ie)) ;
        % determine the overall thermal conductivity in this element
        lambda(ie) = lambda_s*(1 - phi) + lambda_i*phi_i(ie) + lambda_w*(phi - phi_i(ie)) ;
        % local heat storage matrix
        M = C_p(ie)*shapeshape_tri(coord);      
        % local heat transfer matrix 
        S = dshapedshape_tri(coord,[lambda(ie), 0.0; 0.0, lambda(ie)]);
        % local phase change matrix
        P = rho_w*L_I*alpha*shapeshape_tri(coord) ;
        % local advection matrix
        % Adv  = shapedshape_tri( coord, [vel;0.0]);% 
        % local dispersion/diffusion matrix
        % Disp = dshapedshape_tri(coord, [Dt, 0.0; 0.0, 0.01*Dt]);
        % add advection and dispersion matrix t
        % S    = Adv + Disp;  
        % assemble to the LHS
        l_lhs= ((1.0/dt)*M + theta * S + theta * P); 
        LHS(sctr,sctr) = LHS(sctr,sctr) + l_lhs; 
        %phi_i_global(ie) =  phi_i(ie);
        %T_ele_global(ie) = T_ele(ie);
        % assemble to the RHS
        l_rhs= ( (1.0/dt)*M - ((1-theta) * S) +  (1 - theta) * P)* u_pre(sctr); %+ (1 - theta) * P * ones(3,1);  
        RHS(sctr) = RHS(sctr) + l_rhs; 
    end


    
    % check out the sparsity of LHS
    % spy(LHS);
    % imposing boundary conditions---------------------
    % loop over all boundary condition nodes
    for ib = 1 : nbcn
        idx = index_bcn(ib);
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
    u_cur = bicgstab(LHS,RHS,1e-6,max_lin_sol_iter) ; % solve the linear system
    % u_cur = LHS \ RHS;
    error_list = u_cur - u_comparion ;
    error = max(error_list) ;
    counter = counter + 1 ;
    
  end
    % copy current result to previous 
    % also initial guess of next time step
    u_pre = u_cur; 
    phi_i_old = phi_i;
    % save current result
    u_record(:,ti+1) = u_cur;  
    phi_i_record(:,ti+1) = phi_i;
    T_ele_record(:,ti+1) = T_ele;
    C_p_record(:,ti+1) = C_p;
    lambda_record(:,ti+1) = lambda;
    str2 = ['iteration times ', num2str(counter)];
    disp(str2);
    
    
    
    % hold on;
end

% Postprocessing-------------------------------
% 
% hold on;

% take the right boundary and plot their values. 
index_bc_right = find(nodes(:,1) == 1.0);
bc_right_nodes = nodes(index_bc_right,:);
bc_right_values= u_record(index_bc_right, ti);

% scatter(nodes(:,1),nodes(:,2),5,u_record(:,10));











