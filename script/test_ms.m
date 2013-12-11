
[H, W, ~] = size(img);
canvas = zeros(H*W, 3);

mn = mean(y(1:3:end));
sd = std(y(1:3:end));
canvas(idx, 1) = (y(1:3:end)-mn)/sd+1;

mn = mean(y(2:3:end));
sd = std(y(2:3:end));

canvas(idx, 2) = (y(2:3:end)-mn)/sd+1;

mn = mean(y(3:3:end));
sd = std(y(3:3:end));
canvas(idx, 3) = (y(3:3:end)-mn)/sd+1;
canvas = reshape(canvas, [H, W, 3]);

subplot(1,2,1);
imshow(img);
subplot(1,2,2);
imshow(canvas);

% [seg] = meanShiftPixCluster(canvas, 1, 0.2, 0.25, false);
% 
% subplot(2,2,1);
% imshow(img);
% 
% subplot(2,2,2);
% imagesc(seg(:,:,1)); axis equal;
% 
% subplot(2,2,3);
% imagesc(seg(:,:,2)); axis equal;
% 
% subplot(2,2,4);
% imagesc(seg(:,:,3)); axis equal;