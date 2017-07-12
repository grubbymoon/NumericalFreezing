% shape dshape product for line element
% input parameters
% coordinates: 2 * 2 matrix, x and y coordinates of connected nodes
% coeff      : 2 * 1 vector, coefficient vector in x and y direction
% output parameters
% dNTN       : 2 * 2 matrix, values of dN'*N 

function NTdN = shapedshape_line(coeff)

% should be 2*1 
%[nrows,ncols] = size(coordinates);

% get area
%I = ones(nrows,1);
%A = 0.5 * det([I coordinates]);

% this is dshapedshape
% dNTdN = A * [-1, -1; 
%              1, 1]; 

%dN = A * [ -1 , 1 ];  

% dNTdN = A * 100* [-1, 1; 
%              -1, 1]; 

%dNTdN = dN' * coeff * dN; 
NTdN = coeff/2*[1 -1;-1,1];