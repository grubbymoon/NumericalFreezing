 % % % % % % % % % % % % % % % % % % % % % % % % 
% freezing equation solver without mechanics 1D (Newton method for
% nonlinearlarity) compared with comsol result with latent heat
% no porosity and all the water is in solid phase initially sigmoid
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
rho_w = 997 ; % kg*m3-1
% define heat capacity of ice 
C_i = 2108 ;
% define thermal conductivity of ice
lambda_i = 2.14 ;
% define density of ice
rho_i = 917 ; 
% define real density of ice (since we neglect volume change and mechanics
% of ice) we set the density of ice and water same)
rho_IR = 917 ;
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
phi = 0.05 ;
% define the alpha 
% alpha = -1e-3;
% alpha = -10 ;
w = 20;
% define the overall heat capacity in one element 
% define the overall thermal conductivity in one element

% boundary temperature
%T_B_amb = -10;     % Celsius degree
%T_A_in  = 10;     % Celsius degree
% end of parameter input ----------------------

% Preprocessing -------------------------------
% reading the nodes
nodes      = dlmread('nodes_coordinates_line.txt');
nodes      = nodes(:,:) ;
% reading the element connectivity
elements   = dlmread('elem_connectivity_line.txt');
% convert to matlab index start from one
elements   = elements + 1; 
% now the spacial dimension and element type
[nn,sdim]  = size(nodes);
[ne,dummy] = size(elements);
% plot the mesh to check

% boundary conditions
index1 = find(nodes(:,1) == 0); % get the index of nodes that the x value equals 0
%index2 = find(nodes(:,2) > 0.075); % get the index of nodes that the y value bigger than 0.075
%index3 = find(nodes(:,2) < 0.125); % get the index of nodes that the y value smaller than 0.125
% index of the three regions above 
% index_bcn  = intersect(intersect(index1, index2), index3); 
index_bcn = index1;
% assign the x,y value of the boundary into 'bc_middle_nodes'
bc_middel_nodes = nodes(index_bcn,:);
% assign the number of row and column of 'bc_middle_nodes' to 'nbcn' and 'dummy'
[nbcn,dummy] = size( bc_middel_nodes );
bc_nodes_value = -6 * ones(nbcn,1);
% bc_ambient_nodes are the rest part
% index = find(bc_all_nodes ~= bc_middel_nodes);
% bc_ambient_nodes = bc_all_nodes(index,:);

% initial conditions (note! T is for every node, other parameters are for elements)
T_ini = ones(nn,1); 
T_ini = 4 * T_ini; 
phi_i_initial = 0*ones(ne,1); 
C_p_initial = C_i*ones(ne,1);
lambda_initial = lambda_i*ones(ne,1);

%phi_w_inital = 0.3*ones(ne,1);


% Processing ----------------------------------
% numerical control
max_lin_sol_iter = 1000; 
% time control
% time_steps = [0:10:400, 500:100:10000, 11000:1000:86400]';
% time_steps = [0:10:400, 500:100:10000]';
% time_steps = [0:1:1200]';
% time_steps = [0:864:8640]'; % time interval is 864s for every step
% time_steps = [0:10:100]';
 time_steps = [0:900:86400]';
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
u_cur(1) = -6;

% initialize parameters
T_ele = sparse(ne,1 );
T_ele_j = sparse(ne,1 );
phi_i = sparse(ne,1 );
phi_i_j = sparse(ne,1 );
C_p = sparse(ne,1 );
C_p_j = sparse(ne,1 );
lambda = sparse(ne,1);
lambda_j = sparse(ne,1);
%phi_w = sparse(ne,1 );
phi_i_j = phi_i_initial;
phi_i_old_j = phi_i_initial;
phi_i_old = phi_i_initial;
phi_i = phi_i_initial;

C_p = C_p_initial;
C_p_j = C_p_initial;
lambda = lambda_initial;
lambda_j = lambda_initial;
% phi_w = phi_w_inital;
% tol = 1e-4 ; % set the tolerance of picard
counter = 0 ; % count the times of iteration
maxiter = 100 ;
% set the increment
h = 1e-4;
tol = 1e-4;

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
% change the initial node temperature to initial element temperature
for ts = 1 : ne 
    noele = elements(ts,:);
    T_ele_record(ts,1) = mean(u_cur(noele));
end
% get the x,y coordinate for every element
% coorele = zeros(ne,2)
% for cs = 1 : ne
  %  coorele(cs,1) = (nodes(elements(cs,1),1) + nodes(elements(cs,2),1) + nodes(elements(cs,3),1))/3;
  %  coorele(cs,2) = (nodes(elements(cs,1),2) + nodes(elements(cs,2),2) + nodes(elements(cs,3),2))/3;
% end 



    

% time step loop 
for ti = 1 : steps
    % cleaning of LHS and RHS
    LHS    = sparse(nn,nn);
    RHS    = sparse(nn,1 );
    phi_i = sparse(ne,1);
    T_ele = sparse(ne,1);
    C_p = sparse(ne,1 );
    lambda = sparse(ne,1);
    LHS_j    = sparse(nn,nn);
    RHS_j    = sparse(nn,1 );
    phi_i_j = sparse(ne,1);
    T_ele_j = sparse(ne,1);
    C_p_j = sparse(ne,1 );
    lambda_j = sparse(ne,1);
    % cleaning of phi
    str = ['Time step ', num2str(ti), ' at time ', num2str(time_steps(ti+1))];
    disp(str);
    dt = time_steps(ti+1) - time_steps(ti); 
    res = 1 ; % set the initial residual
    counter = 0 ;
    
   while (res > tol) && (counter < maxiter)
    
    LHS    = sparse(nn,nn);
    RHS    = sparse(nn,1 );
    LHS_j    = sparse(nn,nn);
    RHS_j    = sparse(nn,1 );
    % u_old = u_cur;
    counter = counter + 1 ;
    u_norm = norm(u_cur);
    
    % loop over all the elements, 
    for ie = 1 : ne % 'ne' is the number of the total elements 
        % get the coordinates of connecting nodes
        sctr  = elements(ie,:) ; % get the nodes of the local element 
        coord = nodes(sctr,:) ; % get the coordinates of the nodes 
        % initial guess of the temperature of the element based on nodes
        % u_comparison = u_cur;
        T_ele(ie) = mean(u_cur(sctr)) ;
      
        T_ele_j(ie) = mean(u_cur(sctr)) + u_norm*h ;
        % evaluate the freezing value
        % rhophrate_i = freezing(T, phi_i(ie));
        switch 2
            case 1
        if phi_i_old(ie) <= 0 && T_ele(ie) > T_m ;
            rhophrate_i = 0; % when there is no ice, no more melting
        elseif phi - phi_i_old(ie) <= 0 && T_ele(ie) < T_m ;
            rhophrate_i = 0; % when there is no water, no more freezing
        else 
            rhophrate_i = alpha*(T_ele(ie) - T_m) ;
        end
        
        if phi_i_old_j(ie) <= 0 && T_ele_j(ie) > T_m ;
            rhophrate_i_j = 0; % when there is no ice, no more melting
        elseif phi - phi_i_old_j(ie) <= 0 && T_ele_j(ie) < T_m ;
            rhophrate_i_j = 0; % when there is no water, no more freezing
        else 
            rhophrate_i_j = alpha*(T_ele_j(ie) - T_m) ;
        end
        % evaluate the change of porosity
        phi_i(ie) = phi_i_old(ie) + dt*rhophrate_i/rho_IR;
        phi_i_j(ie) = phi_i_old_j(ie) + dt*rhophrate_i_j/rho_IR;
        % make a limitation of the ice volume fraction
        phi_i(ie) = min(phi_i(ie),phi);
        phi_i(ie) = max(phi_i(ie),0);
        phi_i_j(ie) = min(phi_i_j(ie),phi);
        phi_i_j(ie) = max(phi_i_j(ie),0);
        % determine the overall heat capacity in this element
        C_p(ie) = rho_s*C_s*(1 - phi) + rho_i*C_i*phi_i(ie) + rho_w*C_w*(phi - phi_i(ie)) ;
        C_p_j(ie) = rho_s*C_s*(1 - phi) + rho_i*C_i*phi_i_j(ie) + rho_w*C_w*(phi - phi_i_j(ie)) ;
        % determine the overall thermal conductivity in this element
        lambda(ie) = lambda_s*(1 - phi) + lambda_i*phi_i(ie) + lambda_w*(phi - phi_i(ie)) ;
        lambda_j(ie) = lambda_s*(1 - phi) + lambda_i*phi_i_j(ie) + lambda_w*(phi - phi_i_j(ie)) ;
        % local heat storage matrix
        M = C_p(ie)*shapeshape_line_lump(coord);   
        M_j = C_p_j(ie)*shapeshape_line_lump(coord); 
        % local heat transfer matrix 
        S = dshapedshape_line(coord,lambda(ie));
        S_j = dshapedshape_line(coord,lambda_j(ie));
        % local phase change matrix
        P = L_I*rhophrate_i*shapeshape_line_consistent(coord) ;
        P_j = L_I*rhophrate_i_j*shapeshape_line_consistent(coord) ;
        % P = L_I*alpha*shapeshape_tri(coord) ;
        % local advection matrix
        % Adv  = shapedshape_tri( coord, [vel;0.0]);% 
        % local dispersion/diffusion matrix
        % Disp = dshapedshape_tri(coord, [Dt, 0.0; 0.0, 0.01*Dt]);
        % add advection and dispersion matrix t
        % S    = Adv + Disp;  
        % assemble to the LHS
        l_lhs= ((1.0/dt)*M + theta * S - theta * P); 
        l_lhs_j= ((1.0/dt)*M_j + theta * S_j - theta * P_j); 
        LHS(sctr,sctr) = LHS(sctr,sctr) + l_lhs; 
        LHS_j(sctr,sctr) = LHS_j(sctr,sctr) + l_lhs_j; 
        % phi_i_global(ie) =  phi_i(ie);
        % T_ele_global(ie) = T_ele(ie);
        % assemble to the RHS
        l_rhs= ( (1.0/dt)*M - ((1-theta) * S) - (1-theta) * P)* u_pre(sctr);  
        l_rhs_j= ( (1.0/dt)*M_j - ((1-theta) * S_j) - (1-theta) * P_j)* u_pre(sctr); 
        RHS(sctr) = RHS(sctr) + l_rhs; 
        RHS_j(sctr) = RHS_j(sctr) + l_rhs_j; 
            case 2
         logistic = 1/(1 + exp( -w*T_ele(ie)));
         logistic_j = 1/(1 + exp( -w*T_ele_j(ie)));
         phi_i(ie) = phi*(1 - logistic);
         phi_i_j(ie) = phi*(1 - logistic_j);
         C_p(ie) = rho_s*C_s*(1 - phi) + rho_i*C_i*phi_i(ie) + rho_w*C_w*(phi - phi_i(ie)) ;
         C_p_j(ie) = rho_s*C_s*(1 - phi) + rho_i*C_i*phi_i_j(ie) + rho_w*C_w*(phi - phi_i_j(ie)) ;
         lambda(ie) = lambda_s*(1 - phi) + lambda_i*phi_i(ie) + lambda_w*(phi - phi_i(ie)) ;
         lambda_j(ie) = lambda_s*(1 - phi) + lambda_i*phi_i_j(ie) + lambda_w*(phi - phi_i_j(ie)) ;
         sigmoid_dot = - phi*w*logistic*(1 - logistic) ;
         sigmoid_dot_j = - phi*w*logistic_j*(1 - logistic_j) ;
         M = (C_p(ie) - rho_i*sigmoid_dot*L_I)*shapeshape_line_lump(coord);
         M_j = (C_p_j(ie) - rho_i*sigmoid_dot_j*L_I)*shapeshape_line_lump(coord);
         S = dshapedshape_line(coord,lambda(ie));
         S_j = dshapedshape_line(coord,lambda_j(ie));
         l_lhs = ((1.0/dt)*M + theta * S ); 
         l_lhs_j = ((1.0/dt)*M_j + theta * S_j ); 
         LHS(sctr,sctr) = LHS(sctr,sctr) + l_lhs; 
         LHS_j(sctr,sctr) = LHS_j(sctr,sctr) + l_lhs_j; 
         l_rhs= ( (1.0/dt)*M - ((1-theta) * S))* u_pre(sctr);  
         l_rhs_j= ( (1.0/dt)*M_j - ((1-theta) * S_j))* u_pre(sctr); 
         RHS(sctr) = RHS(sctr) + l_rhs; 
         RHS_j(sctr) = RHS_j(sctr) + l_rhs_j;
            end

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
    
    for ib = 1 : nbcn
        idx = index_bcn(ib);
        % b(i) -= A(i,ii)*u_bar
        RHS_j = RHS_j - LHS_j(:,idx) * bc_nodes_value(ib);
        % A(ii,ii) -> xii
        xii = LHS_j(idx,idx);
        % A(ii,j ) =  0
        LHS_j(idx,:) = 0.0;
        % A(i ,ii) = 0
        LHS_j(:,idx) = 0.0;
        % b(ii) = xii * u_bar
        RHS_j(idx) = xii * bc_nodes_value(ib);
        % A(ii,ii) <- xii
        LHS_j(idx,idx) = xii; 
    end
    
    residual = RHS - LHS*u_cur;
    res = norm(residual);
    disp(res) ;
    if(res < tol)
        break
    end  
    J = calc_Jacobi(u_cur, LHS_j, RHS_j, LHS, RHS, residual, h, u_norm);
    %impose boundary condition on Jacobian Matrix
   for ib = 1 : nbcn
            idx = index_bcn(ib);
                 
            j_i_i = J(idx,idx);
          
            J(idx,:) = 0.0;
          
            J(:,idx) = 0.0;
          
            residual(idx) = 0;
           
           J(idx,idx) = 1;
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
    % u_cur = bicgstab(LHS,RHS,1e-6,max_lin_sol_iter) ; % solve the linear system
    % u_cur = LHS \ RHS;

    

    dT = J \ -residual;
    % dT(101) = bc_nodes_value ;
    u_cur = u_cur + dT;
    % disp(dT);
  end
    % copy current result to previous 
    % also initial guess of next time step
    u_pre = u_cur; 
  for ts = 1 : ne % change the initial node temperature to initial element temperature
    noele = elements(ts,:);
    T_ele_record(ts,ti+1) = mean(u_cur(noele));
  end
    phi_i_old = phi_i;
    phi_i_old_j = phi_i_j;
    % save current result
    u_record(:,ti+1) = u_cur;  
    phi_i_record(:,ti+1) = phi_i;
    %T_ele_record(:,ti+1) = T_ele;
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
% plot the temperature versus distance value for the first three time steps
%plot(nodes(:,1),u_record(:,2),nodes(:,1),u_record(:,3),nodes(:,1),u_record(:,4))
%legend('10s','20s','30s')
%xlabel('distance')
%ylabel('temperature')
%title('alpha = -1e-2')



% scatter(nodes(:,1),nodes(:,2),5,u_record(:,10));





% plot(nodes(:,1),u_record(:,1),nodes(:,1),u_record(:,16),nodes(:,1),u_record(:,31),nodes(:,1),u_record(:,46),nodes(:,1),u_record(:,61),nodes(:,1),u_record(:,121),nodes(:,1),u_record(:,181),nodes(:,1),u_record(:,241),nodes(:,1),u_record(:,301),nodes(:,1),u_record(:,361),nodes(:,1),u_record(:,421),nodes(:,1),u_record(:,481),nodes(:,1),u_record(:,541),nodes(:,1),u_record(:,601),nodes(:,1),u_record(:,661),nodes(:,1),u_record(:,721),nodes(:,1),u_record(:,781),nodes(:,1),u_record(:,841),nodes(:,1),u_record(:,901),nodes(:,1),u_record(:,961),nodes(:,1),u_record(:,1021),nodes(:,1),u_record(:,1081),nodes(:,1),u_record(:,1141),nodes(:,1),u_record(:,1201))





