function [ Matrix ] = Vector_Combine( A,B )
%A B is n*1 vector
%Output Matrix is a 2n*1 vector
Comp={A;B};
Matrix=cell2mat(Comp);


end

