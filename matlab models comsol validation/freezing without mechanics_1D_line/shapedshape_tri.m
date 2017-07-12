% dshape shape product for triangle element
% input parameters
% coordinates: 3 * 2 matrix, x and y coordinates of connected nodes
% coeff      : 2 * 1 vector, coefficient vector in x and y direction
% output parameters
% dNTN       : 3 * 3 matrix, values of dN'*N 

function NTdN = shapedshape_tri(coordinates, coeff)

% should be 3*2 
[nrows,ncols] = size(coordinates);

% coefficient vector
cx = coeff(1); cy = coeff(2);
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


NTdN = [cx*b1+cy*c1, cx*b2+cy*c2, cx*b3+cy*c3; 
        cx*b1+cy*c1, cx*b2+cy*c2, cx*b3+cy*c3;
        cx*b1+cy*c1, cx*b2+cy*c2, cx*b3+cy*c3  ];

NTdN = NTdN * A / 3.0;


