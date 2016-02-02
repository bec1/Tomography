% % % 
%  A quick script that orbits the camera around the 3D figure and prints 
%  frames for an animation
% % % 

n      = 60;               % number of frames 
dtheta = 6;                % azimuthal angle displacement
dphi   = 0;                % altitude angle displacement
animdir= 'frames140ms/';   % directory for the animation frames

for i=1:n
    camorbit(dtheta,dphi)  
    fname=strcat(animdir,num2str(i)); % filename
    drawnow 
    print('-djpeg','-r100',fname) 
end