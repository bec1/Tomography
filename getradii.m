function [rads,xfits] = getradii(raw)

s= size(raw);


for i=1:s(1)
    data = raw(i,:);
    [rads(i),xfits(i).fit] = fitrad(data);
end
end
    
function [rad,xft] = fitrad(data)
if(sum(data) >1)

    M=13/9; %Magnification
    points=M *(1:length(data))'; %x pixel positions in microns


    xopts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    xopts.Display = 'Off';

            xft = fittype( 'A*r*sin(real(acos((x-x0)/r)))+c', 'independent', 'x', 'dependent', 'y' );
            xopts.StartPoint = [1 0 120 400];



    % Fit model to data.
    xfitresult = fit(points, data', xft, xopts );
    rad = xfitresult.r;
else
    rad = 1e3;
    xft = '';
end
end



