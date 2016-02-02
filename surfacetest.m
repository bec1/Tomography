% The sample isosurface plotter

% x,y,z are coordinate matrices. v is your data. all are the same size and
% 3d
[x,y,z,v] = flow;

p = patch(isosurface(x,y,z,v,-3));
isonormals(x,y,z,v,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1,1,1])
view(3); axis tight
camlight 
lighting gouraud