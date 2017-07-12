% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Function: converting MATLAB file to VTK format.                         %
%                                                                         %
% Authors: Reza Zolfaghari and Dmitri Naumov                              %
%                                                                         %
% Email  : Reza.Zolfaghari@ufz.de                                         %
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%matlab2vtk converter
%  xs and ys are coordinates.
%  ts are triangles.
%  fn is filename.
%  cs contains values in points.
% I_tot is the total number of mobile compounds. for ploting immobile species
% the correct number of species should be passed on.

function[] = matlab2vtk(xs,ys,ts,fn,cs,I_tot,nnodes,ne,legen) 
fd = fopen(fn, 'w');
fprintf(fd, '# vtk DataFile Version 3.1\n');
fprintf(fd, 'concentrations.\n');
fprintf(fd, 'ASCII\n');
fprintf(fd, 'DATASET UNSTRUCTURED_GRID\n');

%  Points.
fprintf(fd, '\n');
fprintf(fd, 'POINTS %d FLOAT\n', nnodes);
for i=1:nnodes
    fprintf(fd, '%g %g 0\n', xs(i),ys(i));
end


%  Triangles.
fprintf(fd, '\n');
fprintf(fd, 'CELLS %d %d\n', ne, 4*ne);
for i=1:ne
    fprintf(fd, '3 %d %d %d\n', ts(i,1), ts(i,2), ts(i,3));
end

fprintf(fd, '\n');
fprintf(fd, 'CELL_TYPES %d\n', ne);
for i=1:ne
    fprintf(fd, '5 ');
    
    if (mod(i,10) == 0)
        fprintf(fd, '\n');

    end
    
 end
%  Data.

fprintf(fd, '\n');
fprintf(fd, 'POINT_DATA %d\n', nnodes);


for j = 1:I_tot
    fprintf(fd, 'SCALARS %s FLOAT\n', legen(1,j));
    fprintf(fd, 'LOOKUP_TABLE default\n');


    for i=1:nnodes
       fprintf(fd, '%6.8e\n', cs(j,i));
    end
end

fclose(fd);

