% shape function for line element (lumped element formulation)
% input parameters
% coordinates: 2 * 2 matrix, x and y coordinates of connected nodes
% output parameters
% NTN        : 2 * 2 matrix, values of N'*N 

function NTN = shapeshape_line_lump(coordinates, L)

% should be 2*1 
%[nrows,ncols] = size(coordinates);

%I = ones(nrows,1);

%A = 0.5 * det([I coordinates]);

%NTN = A * [2/3, 1/3;
 %    1/3, 2/3]; 
 

NTN = L/2*[1 0;0 1];