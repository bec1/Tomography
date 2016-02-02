[x,y,z,v] = flow;

p = patch(isosurface(v,0));
isonormals(v,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1,1,1])
view(3); axis tight
camlight 
lighting gouraud