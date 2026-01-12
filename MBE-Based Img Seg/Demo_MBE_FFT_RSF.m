clear all;close all;
c0 = 3;
imgID = 1; % 1,2,3,4,5  % choose one of the five test images

Img = imread([num2str(imgID),'.bmp']);
Img = double(Img(:,:,1));

switch imgID
    case 1
        iterNum = 3000;
        lambda1 = 1.0;
        lambda2 = 2.0;
        nu = 0.003*255*255;% coefficient of the length term
        initialLSF = ones(size(Img(:,:,1))).*c0;
        initialLSF(20:70,30:90) = -c0;
    case 2
        iterNum = 800;
        lambda1 = 0;
        lambda2 = 0;
        nu = 0.003*255*255;% coefficient of the length term
        initialLSF = ones(size(Img(:,:,1))).*c0;
        %         initialLSF(26:32,28:34) = -c0;
    case 3
        iterNum = 1500;
        lambda1 = 1.0;
        lambda2 = 1.0;
        nu = 0.01*255*255;% coefficient of the length term
        initialLSF = ones(size(Img(:,:,1))).*c0;
        initialLSF(15:78,32:95) = -c0;
    case 4
        iterNum =150;
        lambda1 = 1.0;
        lambda2 = 1.0;
        nu = 0.001*255*255;% coefficient of the length term
        initialLSF = ones(size(Img(:,:,1))).*c0;
        initialLSF(53:77,46:70) = -c0;
    case 111
        iterNum =220;
        lambda1 = 1.0;
        lambda2 = 1.0;
        nu = 0.001*255*255;% coefficient of the length term
        initialLSF = ones(size(Img(:,:,1))).*c0;
        initialLSF(47:60,86:99) = -c0;
          case 6
        iterNum =1200;
        lambda1 = 0.6;
       lambda2 = 0.4;
       nu = 400;% coefficient of the length term
        initialLSF = ones(size(Img(:,:,1))).*c0;
        initialLSF(50:80,30:60) = -c0;
end

u = initialLSF;
% [nrow,ncol] =size(Img);
% ic=nrow/2;
% jc=ncol/2;
% r=20;
% phi_0 = sdf2circle(nrow,ncol,ic,jc,r);
% [X,Y] = meshgrid(1:ncol, 1:nrow);
% phi = sqrt((X-jc).^2+(Y-ic).^2)-r;
% u = phi;

size(u)

figure;imagesc(Img, [0, 255]);colormap(jet);hold on;axis off,axis equal
title('Initial contour');
[c,h] = contour(u,[0 0],'r');
pause(0.1);

figure,mesh(u);
title('Initial \phi');

timestep = .001;% time step
mu = 1;% coefficient of the level set (distance) regularization term P(\phi)

epsilon = 1.0;% the papramater in the definition of smoothed Dirac function
sigma = 3;    % scale parameter in Gaussian kernel
% Note: A larger scale parameter sigma, such as sigma=10, would make the LBF algorithm more robust
%       to initialization, but the segmentation result may not be as accurate as using
%       a small sigma when there is severe intensity inhomogeneity in the image. If the intensity
%       inhomogeneity is not severe, a relatively larger sigma can be used to increase the robustness of the LBF
%       algorithm.
K=fspecial('gaussian',round(2*sigma)*2+1,sigma);     % the Gaussian kernel
I = Img;
KI=conv2(Img,K,'same');     % compute the convolution of the image with the Gaussian kernel outside the iteration
% See Section IV-A in the above IEEE TIP paper for implementation.

KONE=conv2(ones(size(Img)),K,'same');  % compute the convolution of Gassian kernel and constant 1 outside the iteration
% See Section IV-A in the above IEEE TIP paper for implementation.

delta = 5;
% start level set evolution
for n=1:iterNum
    u = Reg_RSF_MBE_Fourier(u,I,K,KI,KONE,nu,delta,timestep,1,mu,lambda1,lambda2,epsilon,n);
    if mod(n,50)==0
        pause(0.1);
       % figure,
        % mesh(u)
        % iterNum=[num2str(n), ' iterations'];
        % title(iterNum);

        [gx,gy]=gradient(u);
        gxy = sqrt(gx.^2 + gy.^2);
       % figure,
       % imagesc(gxy);colormap(jet);hold on;axis off,axis equal
       % [c,h] = contour(u,[0 0],'r'); colorbar
       % iterNum=[num2str(n), ' iterations, |\nabla \phi|'];
       % title(iterNum);

        hold off;
    end
end

figure,
imagesc(Img, [0, 255]);colormap(gray);hold on;axis off,axis equal
[c,h] = contour(u,[0 0],'r');
totalIterNum=[num2str(n), ' iterations'];
title(['Final contour, ', totalIterNum]);

figure;
mesh(u);
title('Final level set function');

figure,
[gx,gy]=gradient(u);
gxy = sqrt(gx.^2 + gy.^2);
imagesc(gxy);colormap(jet);hold on;axis off,axis equal
[c,h] = contour(u,[0 0],'r');
colorbar
title('|\nabla \phi|')
save(['MBE_FFT_RSF-' num2str(imgID)]);
