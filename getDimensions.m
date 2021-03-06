function [widths,heights] = getDimensions(cloud,surface_level,zpositions)
    M=(16/15); %Magnification
    img = cloud(:,:,1);
    points=M *(1:size(img,1))'; %x pixel positions in microns
    w=10; %width of slice
    w=floor(w/2);
    for i=1:size(cloud,3)
        img = cloud(:,:,i);
        %img = im2bw(img,surface_level);
        r0=floor(centerfinder(img));
        yslice=mean(img(:,(r0(1)-w):(r0(1)+w)),2);
        xslice=mean(img((r0(2)-w):(r0(2)+w),:),1);

        % Set up fittype and options.
        ft = fittype( 'a*erf((x-x0)/(sqrt(2)*w0))+a*erf(-(x-x1)/(sqrt(2)*w1))+c', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [1 0 3 3 20 100];

        % Fit model to data.
        xfitresult = fit(points, xslice', ft, opts );

        % Fit model to data.
        yfitresult = fit(points, yslice, ft, opts );
        
        confint_y = confint(yfitresult);
        confint_x = confint(xfitresult);
        
        errheights(:,i) = sqrt((confint_y(1,end-1) - yfitresult.x0).^2 +(confint_y(1,end) - yfitresult.x1).^2) ;
        errwidths(:,i) = sqrt((confint_x(1,end-1) - xfitresult.x0).^2 +(confint_x(1,end) - xfitresult.x1).^2) ;
        heights(:,i) = yfitresult.x1-yfitresult.x0;
        widths(:,i) = xfitresult.x1-xfitresult.x0;

    end
    
    ft = fittype( 'poly1' );% Fit model to data.
    [fitresult, gof] = fit( zpositions', widths', ft );

    figure(5)
    % Plot widths
    subplot(1,2,1)
    errorbar(zpositions,widths,errwidths,'.','MarkerSize',20)
%     hold all 
%     plot(fitresult)
    title('widths')
    xlabel('z (\mum)')
    ylabel('x (\mum)')

%     ft = fittype( 'poly1' );% Fit model to data.
%     [fitresult, gof] = fit( zpositions', heights', ft );
    
    % Plot heights
    subplot(1,2,2)
    errorbar(zpositions,heights,errheights,'.','MarkerSize',20)
%     hold all
%     plot(fitresult)
    title('heights')
    xlabel('z (\mum)')
    ylabel('y (\mum)')

    
    figure(6)
    %eccentricity = sqrt(1-(heights.^2)./(widths.^2));
        plot(zpositions,heights./widths,'.','MarkerSize',20)
%     hold all
%     plot(fitresult)
    title('heights')
    xlabel('z (\mum)')
    ylabel('y/x')
    
end

function r=centerfinder(img)
    x=(1:size(img,1))';
    hsum=sum(img,1)'; 
    vsum=sum(img,2);
    gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
    startPoints = [3 100 100 0.1];
    xfit = fit(x,hsum,gaussEqn,'Start', startPoints);
    yfit = fit(x,vsum,gaussEqn,'Start', startPoints);
    r=[xfit.b,yfit.b];

end