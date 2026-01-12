clear; close all
% name = 'DRLSE_flower_spot_0524';MBE_FFT_
name = 'MBE_SAV_RSF-100007';
load([name '.mat'])

[Nx, Ny] = size(Img);
% 
 p00=figure,
 imagesc(Img,[0, 255]);
  axis off;
 axis equal; 
 colormap(gray);
hold on;
  saveas(p00, [name '_init.png'])

p0 = figure,
mesh(initialLSF);
set(gca,'ZDir','reverse');
hold on;
contour(initialLSF, [0,0], 'r', 'LineWidth', 2);
 saveas(p0, [name '_init_mesh.png'])

p_1 = figure,
imagesc(Img,[0, 255]);
 axis off;
axis equal; colormap(gray);
hold on;  contour(initialLSF, [0,0], 'r', 'LineWidth', 2);
  saveas(p_1, [name '_init_contour.png'])
 
 p_1 = figure,
imagesc(Img,[0, 255]);
 axis off;
axis equal; colormap(jet);
hold on;  contour(initialLSF, [0,0], 'r', 'LineWidth', 2);
colorbar;
 saveas(p_1, [name '_init_contour1.png'])

p1 = figure,
imagesc(Img,[0, 255]); axis off;
axis equal; colormap(gray);
hold on;  contour(phi_0, [0,0], 'r', 'LineWidth', 2);
  saveas(p1, [name '_contour.png'])

p2 = figure,
imagesc(gxy); colormap(jet);
axis equal;
axis([1 Ny 1 Nx]);
 axis off;
colorbar;
 saveas(p2, [name '_nablau.png'])

p3 = figure,
mesh(phi_0);hold on;
set(gca,'ZDir','reverse');
contour(phi, [0,0], 'r', 'LineWidth', 2);
 saveas(p3, [name '_u.png'])

p4 = figure,
hold on; contour(phi, [0, 0], 'r', 'LineWidth', 2);
axis equal;
set(gca,'YDir','reverse');
hold on; quiver(gx, gy)