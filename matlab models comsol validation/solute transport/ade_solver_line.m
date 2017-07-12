% % % % % % % % % % % % % % % % % % % % % % % % 
% Advection-Diffusion equation solver with FEM
% 
% Authors: Haibing Shao
% Email  : Haibing.Shao@ufz.de
% 
% % % % % % % % % % % % % % % % % % % % % % % % 

% The governing equation is: 
% dC/dt + v*dC/dx - Dt*d2C/dy2 + R(C) = 0 
clear all;
% parameter input ----------------------------
% define velocity
vel = 1. / 86400. ;   % 1 m/d
% define dispersion term
% Dt  = 2.5e-4 / 86400. ; % 2.5 cm2/d
Dt  = 0.5 * vel ; 
% stoichiometry
% fa=1;fb=1;fc=1;
f_stoi=[1.0;1.0;1.0];
% decay term
k_dec = 0.1 / 86400;  % 0.1 1/d
% boundary concentration
c_B_amb = 2.5e-4;     % mol/l
c_A_in  = 3.3e-4;     % mol/l
% end of parameter input ----------------------

% Preprocessing -------------------------------
% Reading the nodes
% this is 2D triangular mesh
% nodes    = dlmread('nodes_coordinates.txt');
% this is 1D line mesh
nodes      = dlmread('nodes_coordinates_line.txt');
% reading the element connectivity
% this is 2D triangular mesh
% elements   = dlmread('elem_connectivity.txt');
% this is 1D line mesh
elements   = dlmread('elem_connectivity_line.txt');
% convert to matlab index start from one
elements   = elements + 1; 
% now the spacial dimension and element type
[nn,sdim]  = size(nodes);
[ne,dummy] = size(elements);

% boundary conditions
% index1 = find(nodes(:,1) == 0);
index_bcn = find(nodes(:,1) == 0);
% index2 = find(nodes(:,2) > 0.075);
% index3 = find(nodes(:,2) < 0.125);
% index_bcn  = intersect(intersect(index1, index2), index3);
bc_middel_nodes = nodes(index_bcn,:);
[nbcn,dummy] = size( bc_middel_nodes );
bc_nodes_value = 1.0 * ones(nbcn,1);
% bc_ambient_nodes are the rest part
% index = find(bc_all_nodes ~= bc_middel_nodes);
% bc_ambient_nodes = bc_all_nodes(index,:);

% initial conditions
c_ini = ones(nn,1); 
c_ini = 1e-9 * c_ini;  % Why not we give zero? 

% Processing ----------------------------------
% numerical control
max_lin_sol_iter = 100; 
% time control
% time_steps = [0:100:86400]';
time_steps = [0:864:86400]';
[steps,dummy] = size(time_steps);
steps = steps - 1; 
% theta = 1.0;  % implicit
theta = 0.5;  % Crank-Nicolson
% theta = 0.0;  % explicit

% initialize previous and current solution matrix
u_pre  = sparse(nn,1 ); 
u_cur  = sparse(nn,1 ); 
u_pre  = c_ini; 

% storage space of unknowns after each time step
u_record      = zeros(nn, steps+1);
u_record(:,1) = u_pre;  

% time step loop 
for ti = 1 : steps; 
    % cleaning of LHS and RHS
    LHS    = sparse(nn,nn);
    RHS    = sparse(nn,1 );
    str = ['Time step ', num2str(ti), ' at time ', num2str(time_steps(ti+1))];
    disp(str);
    dt = time_steps(ti+1) - time_steps(ti);    
    % loop over all the elements, 
    for ie = 1 : ne
        % get the coordinates of connecting nodes
        sctr  = elements(ie,:);
        coord = nodes(sctr,:) ; 
        % local mass matrix
        M    = shapeshape_line(  coord );          % no coeff
        % local advection matrix
        % Adv  = shapedshape_line( coord )* [vel;0.0];% 
        % Adv  = shapedshape_line( coord )*vel;% 
        % local dispersion/diffusion matrix
        % Disp = dshapedshape_line(coord )*[Dt, 0.0; 0.0, 0.01*Dt];
        Disp = dshapedshape_line(coord , Dt);
        % add advection and dispersion matrix t
        % S    =  Adv + Disp;  
        S    =  Disp;  
        % assemble to the LHS
        l_lhs= ((1.0/dt)*M + theta    * S); 
        LHS(sctr,sctr) = LHS(sctr,sctr) + l_lhs; 
        % assemble to the RHS
        l_rhs= ((1.0/dt)*M - (1-theta)* S) * u_pre(sctr);  
        RHS(sctr) = RHS(sctr) + l_rhs; 
    end

    % guess how does LHS looks like? 
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

    % solve linear equation system. 
    % u_cur = LHS \ RHS;
    u_cur = bicgstab(LHS,RHS,1e-6,max_lin_sol_iter);
    % copy current result to previous
    u_pre = u_cur; 
    % save current result
    u_record(:,ti+1) = u_cur;  
end

% Postprocessing-------------------------------
% plot the value
plot( nodes(:,1), u_record(:,ti+1)); 



