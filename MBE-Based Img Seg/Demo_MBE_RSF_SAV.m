clear all;close all;
c0 = 3;
imgID = 6; 

Img = imread([num2str(imgID),'.bmp']);

Img = double(Img(:,:,1));

[M,N]=size(Img);
switch imgID
    case 1
        pa.iterNum = 3000;
        pa.lambda1 =1;
        pa.lambda2 = 1.95;
        pa.nu = 0.005*255*255;% coefficient of the length term
        initialLSF = ones(size(Img(:,:,1))).*c0;
        initialLSF(14:70,30:92) = -c0;
    case 2
        pa.iterNum = 2800;
        pa.lambda1 = 4.8;
        pa.lambda2 = 4.8;
        pa.nu = 0.010*255*255;% coefficient of the length term
        initialLSF =ones(size(Img(:,:,1))).*c0;
              initialLSF(12:68,18:64) = -c0;
    case 3
        pa.iterNum = 4000;
        pa.lambda1 = 5;
        pa.lambda2 = 4.5;
         pa.nu = 0.015*255*255;
%          pa.nu = 850;
        initialLSF = ones(size(Img(:,:,1))).*c0;
        initialLSF(10:79,40:92) = -c0;
    case 4
        pa.iterNum =1000;
        pa.lambda1 = 5;
        pa.lambda2 = 3;
        pa.nu = 0.01*255*255;% coefficient of the length term
        initialLSF = ones(size(Img(:,:,1))).*c0;
        initialLSF(23:87,46:70) = -c0;
           case 5
        pa.iterNum =1000;
        pa.lambda1 = 5;
        pa.lambda2 = 3;
        pa.nu = 0.01*255*255;% coefficient of the length term
        initialLSF = ones(size(Img(:,:,1))).*c0;
        initialLSF(23:87,46:70) = -c0;
           case 6
        pa.iterNum =3000;
        pa.lambda1 = 0.6;
        pa.lambda2 = 0.4;
        pa.nu = 400;% coefficient of the length term
        initialLSF = ones(size(Img(:,:,1))).*c0;
        initialLSF(50:80,30:60) = -c0;
           case 7
        pa.iterNum =1000;
        pa.lambda1 = 5;
        pa.lambda2 = 3;
        pa.nu = 0.01*255*255;% coefficient of the length term
        initialLSF = ones(size(Img(:,:,1))).*c0;
        initialLSF(23:87,46:70) = -c0;
    case 9
        pa.iterNum =3300;
        pa.lambda1 = 2;
        pa.lambda2 = 4;
        pa.nu = 0.004*255*255;% coefficient of the length term
        initialLSF = ones(size(Img(:,:,1))).*c0;
        initialLSF(27:75,26:79) = -c0;
end

 phi_0 = initialLSF;
 
figure,mesh(phi_0);
title('Initial \phi');
pa.dt = .001;% time step
pa.mu = 1;% coefficient of the level set (distance) regularization term P(\phi)
pa.epsilon = 1.0;% the papramater in the definition of smoothed Dirac function
pa.sigma = 10;    % scale parameter in Gaussian kernel
pa.K=fspecial('gaussian',round(2*pa.sigma)*2+1,pa.sigma);     
I = Img;
% pa.KI=conv2(Img,pa.K,'same'); 
pa.KI = imfilter(Img, pa.K, 'replicate')
pa.KONE=conv2(ones(size(Img)),pa.K,'same');  % compute the convolution of Gassian kernel and constant 1 outside the iteration
pa.delta = 40;
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
% save(['MBE_SAV_RSF-' num2str(imgID)]);