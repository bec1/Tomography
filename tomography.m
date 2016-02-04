function tomography(sideslices,topslices)
%% TOMOGRAPHY loads slice images taken from the top and from the side, and plots them.
%  Usage:   inputs:   
%                   - sideslices = loadfits(side

%% Load data sets
side_data = loadDataset(sideslices);
top_data = loadDataset(topslices);

%% Crop parameters
side_a = 170; % width of square crop
side_crop = [66,304,side_a,side_a];

top_w = 60;  % width of crop (radial)
top_l = 400; % length of crop (axial)
top_r = 0;  % degrees of rotation
top_crop = [225,62,top_w,top_l];

%% Process data
cloud = stackSlices(side_data,side_crop); % crop and stack the slices into a 3dmatrix
cloud_smooth = smoothCloud(cloud); % smooth the cloud for isosurface plotting

profiles = axialProfiles(top_data,top_crop,top_r); % crop and stack the top slices into 1d axial profiles

[zscale,zpositions] = getzs(profiles);

%% Plotting
figure(1)

subplot(2,2,1)
surface_level = 0.2;
p=getIsosurface(cloud_smooth,surface_level,zpositions);
title('Isosurface')

subplot(2,2,2)
hold all;
for i=1:size(profiles,2)
    plot(zscale,profiles(:,i))
end
plot(zscale,sum(profiles,2))
title('all axial profiles')

subplot(2,2,3)
imagesc(cloud(:,:,end/2));
axis image;
title('Sample slice')

subplot(2,2,4)
imagesc(imrotate(imcrop(top_data(end/2).img,top_crop),top_r));
axis image;
title('Sample axial profile')


end

function p=getIsosurface(cloud,surface_level,zpositions)
%% Get and plot the isosurface from a set of slices
    s = size(cloud);
    mag = 13/9;
    x = mag*(1:s(1));
    y = mag*(1:s(2));
    z = zpositions;
    [X,Y,Z] = meshgrid(x,y,z);
    p = patch(isosurface(X,Y,Z,cloud,surface_level));
    isonormals(X,Y,Z,cloud,p)
    p.FaceColor = 'red';
    p.EdgeColor = 'none';
    daspect([1,1,1])
    view(3); axis tight
    camlight 
    lighting gouraud

end

function [points,zpos] = getzs(profiles)
% get the z positions
    s = size(profiles);
    points = (1:s(1))*13/9;
    for i=1:s(2)
        [xData, yData] = prepareCurveData( points', profiles(:,i ));

        % Set up fittype and options.
        ft = fittype( 'gauss1' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.Lower = [-Inf -Inf 0];
        opts.StartPoint = [0.950973163734602 248.444444444444 29.4216810130658];

        % Fit model to data.
        [fitresult, ~] = fit( xData, yData, ft, opts );
        zpos(:,i) = fitresult.b1; % get zpositions from the profiles
    end
end

function cloud_smooth = smoothCloud(cloud)
    s=size(cloud);
    
    for i=1:s(3)
        cloud_smooth(:,:,i) = medfilt2(cloud(:,:,i));
    end
        
end

function profiles = axialProfiles(data,crop,top_rot)
%% Get 1D axial profiles
    for i=1:length(data)
        img = imrotate(imcrop(data(i).img,crop),top_rot);
        profiles(:,i) = mean(img,2);
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