A = [ 1/2, -1/2, 0,   0,    0; 
     -1/2,  3/2, -1,  0,    0; 
       0,   -1, 4/3, -1/3,  0;
       0,    0, -1/3, 2/3, -1/3;
       0,    0,  0,  -1/3, 1/3];
% A = [ 1/2,   0, 0,   0,    0; 
%        0,   3/2, -1,  0,    0; 
%        0,   -1, 4/3, -1/3,  0;
%        0,    0, -1/3, 2/3, -1/3;
%        0,    0,  0,  -1/3, 1/3];

b = [0,0,0,0,0]';
% b = [6,6,0,0,0]';

x = A \ b; 