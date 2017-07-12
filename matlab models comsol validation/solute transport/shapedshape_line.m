% dshape shape product for triangle element
% input parameters
% coordinates: 3 * 2 matrix, x and y coordinates of connected nodes
% coeff      : 2 * 1 vector, coefficient vector in x and y direction
% output parameters
% dNTN       : 3 * 3 matrix, values of dN'*N 

function NTdN = shapedshape_line(coordinates)

% should be 2*1 
[nrows,ncols] = size(coordinates);

I = ones(nrows,1);

A = 0.5 * det([I coordinates]);
% calculate shape function coefficients

% NTdN = [ 1/2, -1/2;
%         -1/2, 1/2] * A ;

NTdN = [-1/2, 1/2;
        -1/2, 1/2] * A ;

