function  sliceplot(M)
s=size(M);
slices=zeros(s(2),s(1));


for i=1:s(3)%6:2:(s(3)-4) for twice height interpolations
    slices=[slices;rot90(M(:,:,i))];
end

imagesc(slices(181:end,:)>0.07)
axis image
end
    