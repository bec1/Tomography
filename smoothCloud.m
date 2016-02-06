function cloud_smooth = smoothCloud(cloud,filtrad)
%% smooth each slice with median filtering to get the edges
    s=size(cloud);
    
    for i=1:s(3)
        cloud_smooth(:,:,i) = medfilt2(cloud(:,:,i),[filtrad,filtrad]); 
    end
        
end