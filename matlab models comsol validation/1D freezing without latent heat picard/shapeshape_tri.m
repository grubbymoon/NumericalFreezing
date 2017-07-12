% shape function for triangle element
% input parameters
% coordinates: 3 * 2 matrix, x and y coordinates of connected nodes
% output parameters
% NTN        : 3 * 3 matrix, values of N'*N 

function NTN = shapeshape_tri(coordinates)

% should be 3*2 
[nrows,ncols] = size(coordinates);

I = ones(nrows,1);

A = 0.5 * det([I coordinates]);

NTN = A * 1./6.* [1    1/2  1/2;
                  1/2  1    1/2;
                  1/2  1/2  1   ]; 


