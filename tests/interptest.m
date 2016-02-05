x(:,1) = [0,1,0,];
x(:,2) = [1,0,1];
x(:,3) = [0,1,0];
subplot(1,3,1)
imagesc(x)
subplot(1,3,2)
y= interp2(x,5);
imagesc(y)
subplot(1,3,3)
y= interp2(x,5,'cubic');
imagesc(y)
