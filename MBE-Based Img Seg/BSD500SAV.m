% BSD500 Image Segmentation with MBE regularizer
% -------------------------------------------------------------------------
% Img = imread('../../TestImages/100007.jpg'); Img = rgb2gray(Img); Img =
% double(Img); Img = imresize(Img, 0.5); initialLSF(70:100, 50:110) = -c0;
% -------------------------------------------------------------------------
% Img = imread('../../TestImages/130014.jpg'); Img = rgb2gray(Img); Img =
% double(Img); Img = imresize(Img, 0.5); initialLSF(112:138, 41:88) = -c0;
% -------------------------------------------------------------------------
% Img = imread('../../TestImages/196062.jpg'); Img = rgb2gray(Img); Img =
% Img(:,[1:248]); Img = double(Img); Img = imresize(Img, 0.5); c0 = 2;
% load('196062_initialLSF') initialLSF = initialLSF(:, 1:124);
% -------------------------------------------------------------------------
% Img = imread('../../TestImages/100080.jpg'); Img = Img(:,:,2); Img =
% double(Img); Img = imresize(Img, 0.5); initialLSF(100:159, 55:109) = -c0;
% -------------------------------------------------------------------------



clear all;close all;
c0 = 2;

Img = imread('196027.jpg'); 

% Img = rgb2gray(Img);
Img = Img(:,:,2);
Img = double(Img); 
Img = imresize(Img, 0.5); 


pa.iterNum = 1200;
pa.lambda1 = 0.40;
pa.lambda2 = 0.60;
pa.nu = 450;%0.001*255*255;% coefficient of the length term
initialLSF = ones(size(Img(:,:,1))).*c0;
% initialLSF(70:100, 50:110) = -c0;
initialLSF(70:180, 65:105) = -c0;
phi_0 = initialLSF;
 
figure,mesh(phi_0);
title('Initial \phi');
pa.dt = .001;% time step
pa.mu = 1;% coefficient of the level set (distance) regularization term P(\phi)
pa.epsilon = 1.0;% the papramater in the definition of smoothed Dirac function
pa.sigma =5;    % scale parameter in Gaussian kernel
pa.K=fspecial('gaussian',round(2*pa.sigma)*2+1,pa.sigma);     
I = Img;
pa.KI=conv2(Img,pa.K,'same');     
pa.KONE=conv2(ones(size(Img)),pa.K,'same');  % compute the convolution of Gassian kernel and constant 1 outside the iteration
pa.delta = 10;
% start level set evolution
figure;imagesc(Img, [0, 255]);colormap(jet);hold on;axis off,axis equal
title('Initial contour');
[c,h] = contour(phi_0,[0 0],'r');
pause(0.1);
for n=1:pa.iterNum
phi_0= SAV_fun(phi_0,Img,pa);
%  if mod(n,100)==0
%         pause(0.1);
%          gx = (phi_0([2:end 1],:)-phi_0([end 1:end-1],:))/2;
%          gy = (phi_0(:,[2:end 1])-phi_0(:,[end 1:end-1]))/2;
%         gxy = sqrt(gx.^2 + gy.^2);
%         figure,
%        imagesc(gxy);colormap(jet);hold on;axis off,axis equal
%        [c,h] = contour(phi_0,[0 0],'r'); colorbar
%        iterNum=[num2str(n), ' iterations, |\nabla \phi|'];
%        title(iterNum);
%         hold off
%  end
end
phi=phi_0;
figure,
imagesc(Img, [0, 255]);colormap(gray);hold on;axis off,axis equal
[c,h] = contour(phi,[0 0],'r');
totalIterNum=[num2str(n), ' iterations'];
title(['Final contour, ', totalIterNum]);
figure;
mesh(phi);
title('Final level set function');
figure,
[gx,gy]=gradient(phi);
gxy = sqrt(gx.^2 + gy.^2);
imagesc(gxy);colormap(jet);hold on;axis off,axis equal
[c,h] = contour(phi,[0 0],'r');
colorbar
title('|\nabla \phi|')
 save(['MBE_SAV_RSF-196027']);