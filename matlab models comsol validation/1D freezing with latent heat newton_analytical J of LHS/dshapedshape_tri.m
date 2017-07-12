% dshape shape product for triangle element
% input parameters
% coordinates: 3 * 2 matrix, x and y coordinates of connected nodes
% coeff      : 2 * 1 vector, coefficient vector in x and y direction
% output parameters
% dNTN       : 3 * 3 matrix, values of dN'*N 

function dNTdN = dshapedshape_tri(coordinates, coeff)

% should be 3*2 
[nrows,ncols] = size(coordinates);

% get area
I = ones(nrows,1);
A = 0.5 * det([I coordinates]);
% calculate shape function coefficients
x1 = coordinates(1,1); 
x2 = coordinates(2,1); 
x3 = coordinates(3,1); 

y1 = coordinates(1,2); 
y2 = coordinates(2,2); 
y3 = coordinates(3,2); 

b1 = (y2 - y3)/2.0/A ; 
b2 = (y3 - y1)/2.0/A ; 
b3 = (y1 - y2)/2.0/A ;

c1 = (x3 - x2)/2.0/A ;
c2 = (x1 - x3)/2.0/A ; 
c3 = (x2 - x1)/2.0/A ; 

% this is dshape
dN = [b1, b2, b3; 
      c1, c2, c3 ];
  
dNTdN = dN' * coeff * dN; 

