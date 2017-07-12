function [ dNTf ] = dshapef_tri( coordinates, coeff )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
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
  dNTf=dN' * coeff;

end

