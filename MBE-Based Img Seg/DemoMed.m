
clear all;close all;
c0 = 2;

Img = imread('01CF3F82AD34B1B7.jpg'); 
Img = Img(187:343, 123:414, 1); 
Img = double(Img); 
% Img = imresize(Img, 0.5); 

pa.iterNum = 1000;
pa.lambda1 = 0.5;
pa.lambda2 = 0.5050;
pa.nu = 200;%0.001*255*255;% coefficient of the length term
initialLSF = ones(size(Img(:,:,1))).*c0;
%  initialLSF(72:129, 141:197) = -c0;
initialLSF(14:144,14:279) = -c0;

pa.dt = 0.0041;% time step
pa.mu = 10;% coefficient of the level set (distance) regularization term P(\phi)
pa.delta = 20;

pa.epsilon = 1;% the papramater in the definition of smoothed Dirac function
pa.sigma =20;    % scale parameter in Gaussian kernel


phi_0 = initialLSF;
% [nrow,ncol] =size(Img);
% ic=nrow/2;
% jc=ncol/2;
% r=50;
% phi_0 = sdf2circle(nrow,ncol,ic,jc,r);
% [X,Y] = meshgrid(1:ncol, 1:nrow);
% phi = sqrt((X-jc).^2+(Y-ic).^2)-r;
% u = phi;


figure;imagesc(Img, [0, 255]);hold on;axis off,axis equal
title('Initial contour');
[c,h] = contour(phi_0,[0 0],'r', 'LineWidth', 1);
pause(0.1);

% a = Img(end:-1:1, end:-1:1);
% b = Img(end:-1:1, :);
% c = Img(:, end:-1:1);
% ba = [b, a];
% dc = [Img, c];
% Img = [ba; dc];
% % cd = [c, Img];
% % ab = [a, b];
% % Img = [Img, c];%[ab; cd];
% 
% ua = u(end:-1:1, end:-1:1);
% ub = u(end:-1:1, :);
% uc = u(:, end:-1:1);
% uba = [ub, ua];
% udc = [u, uc];
% u = [uba; udc];
% % ucd = [uc, u];
% % uab = [ua, ub];
% % u = [u, uc];%[uab; ucd];


pa.K=fspecial('gaussian',round(2*pa.sigma)*2+1,pa.sigma);     % the Gaussian kernel
I = Img;
% KI=conv2(Img,K,'same');     % compute the convolution of the image with the Gaussian kernel outside the iteration
pa.KI = imfilter(Img, pa.K, 'replicate');
% KONE=conv2(ones(size(Img)),K,'same');  % compute the convolution of Gassian kernel and constant 1 outside the iteration
pa.KONE = imfilter(ones(size(Img)), pa.K, 'replicate');
for n=1:pa.iterNum
phi_0= SAV_fun(phi_0,Img,pa);
%  if mod(n,50)==0
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
phi_0=phi_0;
figure,
imagesc(Img, [0, 255]);colormap(gray);hold on;axis off,axis equal
[c,h] = contour(phi_0,[0 0],'r');
totalIterNum=[num2str(n), ' iterations'];
title(['Final contour, ', totalIterNum]);
figure;
mesh(phi_0);
title('Final level set function');
figure,
[gx,gy]=gradient(phi_0);
gxy = sqrt(gx.^2 + gy.^2);
imagesc(gxy);colormap(jet);hold on;axis off,axis equal
[c,h] = contour(phi_0,[0 0],'r');
colorbar
title('|\nabla \phi|')
 save(['MBE_SAV_RSF01CF3F82AD34B1B7' ]);