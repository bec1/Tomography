function cloud_sort = tomography(sideslices,topslices)
%% TOMOGRAPHY loads slice images taken from the top and from the side, and plots them.
%  Usage:   inputs:   
%                   - sideslices = loadfits(side

%% Load data sets
side_data = loadDataset(sideslices);
top_data = loadDataset(topslices);

%% Crop parameters
side_a = 170; % width of square crop
side_crop = [66,304,side_a,side_a];

top_w = 40;  % width of crop (radial)
top_l = 300; % length of crop (axial)
top_r = 2;  % degrees of rotation
top_crop = [245,62,top_w,top_l];

surface_level = 0.5; % Isosurface equidensity level

%% Process data
cloud = stackSlices(side_data,side_crop); % crop and stack the slices into a 3dmatrix
cloud = renormalizeCloud(cloud); % correct for atom number fluctuations - keep on when axial moving
cloud_smooth = smoothCloud(cloud,15); % smooth the cloud for isosurface plotting

profile_imgs = axialProfiles(top_data,top_crop,top_r); % crop and stack the top slices into 1d axial profiles
profiles = squeeze(mean(profile_imgs,2));

[zscale,zpositions] = getzs(profiles); % get the z-positions from the profiles
[zsort,ix] = sort(zpositions); % sort in increasing order
profile_imgs = profile_imgs(:,:,ix); % sort the profiles and clouds
cloud_sort = cloud_smooth(:,:,ix);

getDimensions(cloud_sort,surface_level,zsort); % get the widths and heights

%% Plotting
figure(1)

% Volumetric plot
subplot(2,2,1)
getIsosurface(cloud_sort,surface_level,zsort);
hold all
%plotslice(cloud_smooth,zpositions);
title('Isosurface')

% 1d Axial profiles
subplot(2,2,2)
hold all;
for i=1:size(profiles,2)
    plot(zscale,profiles(:,i))
    % draw a line at the z-positions
    line([zpositions(i) zpositions(i)],[0,1.3*max(profiles(:,i))],'LineWidth',1,'Color','black') 
end
plot(zscale,0.2+sum(profiles,2))
xlim([0.8*min(zpositions), 1.2*max(zpositions)])
title('Axial positions of slices')

% all slices
subplot(2,2,3)
imagesc(reshape(cloud,size(cloud(:,:,1),1),[],1));
axis image;
axis off;
title('All slices so far')

% All the axial profile images
subplot(2,2,4)
imagesc(reshape(profile_imgs,size(profile_imgs(:,:,1),1),[],1));
axis image;
axis off;
title('Axial profiles')


end

function [widths,heights] = getDimensions(cloud,surface_level,zpositions)
    M=(16/15); %Magnification
    img = cloud(:,:,1);
    points=M *(1:size(img,1))'; %x pixel positions in microns
    w=10; %width of slice
    w=floor(w/2);
    for i=1:size(cloud,3)
        img = cloud(:,:,i);
        img = im2bw(img,surface_level);
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

        heights(:,i) = yfitresult.x1-yfitresult.x0;
        widths(:,i) = xfitresult.x1-xfitresult.x0;

    end
    
    ft = fittype( 'poly1' );% Fit model to data.
    [fitresult, gof] = fit( zpositions', widths', ft );

    figure(2)
    % Plot widths
    subplot(1,2,1)
    plot(zpositions,widths,'.','MarkerSize',20)
%     hold all 
%     plot(fitresult)
    title('widths')
    xlabel('z (\mum)')
    ylabel('x (\mum)')

%     ft = fittype( 'poly1' );% Fit model to data.
%     [fitresult, gof] = fit( zpositions', heights', ft );
    
    % Plot heights
    subplot(1,2,2)
    plot(zpositions,heights,'.','MarkerSize',20)
%     hold all
%     plot(fitresult)
    title('heights')
    xlabel('z (\mum)')
    ylabel('y (\mum)')

    
    figure(3)
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

function plotslice(cloud,zpositions)
%% Get and plot the isosurface from a set of slices
    s = size(cloud);
    mag = 1.01;
    x = mag*(1:s(1));
    y = mag*(1:s(2));
    z = zpositions;
    [X,Y,Z] = meshgrid(x,y,z);
%     zq = linspace(min(zpositions),max(zpositions),3*length(zpositions));
%     [Xq,Yq,Zq] = meshgrid(x,y,zq);
%     cloud = interp3(X,Y,Z,cloud,Xq,Yq,Zq);
    slice(X,Y,Z,cloud,100,100,[150,190])
    shading flat
    daspect([1,1,.5])
    view(3); axis tight, box on, grid on, axis vis3d tight
    camlight 
    lighting flat

end

function p=getIsosurface(cloud,surface_level,zpositions)
%% Get and plot the isosurface from a set of slices

    s = size(cloud);
    mag = 1.01;
    x = mag*(1:s(1));
    y = mag*(1:s(2));
    z = zpositions;
    [X,Y,Z] = meshgrid(x,z,y);
    cloud = permute(cloud,[3,2,1]);
    p = patch(isosurface(X,Y,Z,cloud,surface_level));
    isonormals(X,Y,Z,cloud,p)
    p.FaceColor = [222 235 250]/255;
    p.EdgeColor = 'none';
    daspect([1,1,1])
    view(3); axis tight, box on, grid on, axis vis3d tight
    camlight 
    lighting gouraud

end

function [points,zpos] = getzs(profiles)
%% get the z positions
    s = size(profiles);
    points = (1:s(1))*13/9;
    for i=1:s(2)
        [xData, yData] = prepareCurveData( points', profiles(:,i ));
        [~,ix] = max(profiles(:,i));
        startcenter = points(ix);

        % Set up fittype and options.
        ft = fittype( 'gauss1' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.Lower = [-Inf -Inf 0];
        opts.StartPoint = [0.08 startcenter 10];

        % Fit model to data.
        [fitresult, ~] = fit( xData, yData, ft, opts );
        zpos(:,i) = fitresult.b1; % get zpositions from the profiles
    end
end

function cloud_smooth = smoothCloud(cloud,filtrad)
%% smooth each slice with median filtering to get the edges
    s=size(cloud);
    
    for i=1:s(3)
        cloud_smooth(:,:,i) = medfilt2(cloud(:,:,i),[filtrad,filtrad]); 
    end
        
end

function cloudout = renormalizeCloud(cloud)

    for i=1:size(cloud,3)
        cloudout(:,:,i) = 1e4 *cloud(:,:,i)/sum(sum(cloud(:,:,i)));
    end
end

function profiles = axialProfiles(data,crop,top_rot)
%% Get 1D axial profiles
    for i=1:length(data)
        img = imrotate(imcrop(data(i).img,crop),top_rot);
        profiles(:,:,i) = img;
    end

end

function cloud = stackSlices(data,crop)
%% Stack slices and crop them
    for i=1:length(data)
        img = data(i).img;
        cloud(:,:,i) = imcrop(img,crop);
    end
end

function data = loadDataset(img_list)
%% LOAD_DATASET loads the raw images as OD arrays
    % Initialize data struct
    data(1:length(img_list)) = struct('name','','img',[]);
    % Load the images from the filenames
    fprintf('\n');
    for i =1:length(img_list)
        fprintf('.');
        data(i).name = img_list{i};
        data(i).img=loadFitsimage(data(i).name);
    end
    fprintf('\n');
end

function img = loadFitsimage(filename)
 data=fitsread(filename);
    absimg=(data(:,:,2)-data(:,:,3))./(data(:,:,1)-data(:,:,3));

%     % Identify "burned pixels" and make sure they will come out zero.
%     burnedpoints = absimg < 0;
%     absimg(burnedpoints) = 1;
% 
%     % Same thing for points which should be accidentally equal to zero
%     % (withatoms == background) or infinity (withoutatoms == background)
%     zeropoints = absimg == 0;
%     absimg(zeropoints) = 1;
% 
%     infpoints = abs(absimg) == Inf;
%     absimg(infpoints) = 1;
% 
%     nanpoints = isnan(absimg);
%     absimg(nanpoints) = 1;

%replace the pixels with a value of negtive number,0 or inf or nan by the
%average of nearset site.
    ny=size(absimg,1);
    nx=size(absimg,2);
    burnedpoints = absimg <= 0;
    infpoints = abs(absimg) == Inf;
    nanpoints = isnan(absimg);
    Change=or(or(burnedpoints,infpoints),nanpoints);
    NChange=not(Change);
    for i=2:(ny-1)
        for j=2:(nx-1)
            if Change(i,j)
                n=0;
                rp=0;
                if NChange(i-1,j)
                    rp=rp+absimg(i-1,j);
                    n=n+1;
                end
                if NChange(i+1,j)
                    rp=rp+absimg(i+1,j);
                    n=n+1;
                end
                if NChange(i,j-1)
                    rp=rp+absimg(i,j-1);
                    n=n+1;
                end
                if NChange(i,j+1)
                    rp=rp+absimg(i,j+1);
                    n=n+1;
                end
                if (n>0)
                    absimg(i,j)=(rp/n);
                    Change(i,j)=0;
                end
            end
        end
    end
    absimg(Change)=1;
    img = log(absimg);
end