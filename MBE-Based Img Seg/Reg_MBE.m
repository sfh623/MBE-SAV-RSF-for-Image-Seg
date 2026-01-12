function f = Reg_MBE(phi, deta, delta2)

f_lap = phi([1 1:end-1],:) + phi([2:end end],:) + phi(:,[1 1:end-1]) ...
        + phi(:,[2:end end]) - 4*phi;
f_lap2 = f_lap.^2;
f_x = (phi([1 1:end-1],:) - phi([2:end end],:))/2;
f_y = (phi(:,[1 1:end-1]) - phi(:,[2:end end]))/2;
f_g2 = f_x.^2 + f_y.^2;
f2 = (f_g2 - 1).^2;
energy = deta/2*sum(f_lap2(:)) + delta2/4*sum(f2(:));

% [phi_x, phi_y] = gradient(phi);
phi_lap = phi([1 1:end-1],:) + phi([2:end end],:) + phi(:,[1 1:end-1]) ...
    + phi(:,[2:end end]) - 4*phi;
phi_4ord = phi_lap([1 1:end-1],:) + phi_lap([2:end end],:) ...
    + phi_lap(:,[1 1:end-1]) + phi_lap(:,[2:end end]) - 4*phi_lap;
f = -1*deta*phi_4ord;

% % O0 = (abs(phi)>=0.1);
phi_x = (phi([1 1:end-1],:) - phi([2:end end],:))/2;
phi_y = (phi(:,[1 1:end-1]) - phi(:,[2:end end]))/2;
phi_g2 = phi_x.^2 + phi_y.^2;
d_x = (phi_g2-1).*phi_x;
d_y = (phi_g2-1).*phi_y;
dd_x = (d_x([1 1:end-1],:) - d_x([2:end end],:))/2;
dd_y = (d_y(:,[1 1:end-1]) - d_y(:,[2:end end]))/2;

f = f + delta2*(dd_x+dd_y);%.*O0;

phi_mbe = phi + 0.1*f;
f_lap = phi_mbe([1 1:end-1],:) + phi_mbe([2:end end],:) + phi_mbe(:,[1 1:end-1]) ...
        + phi_mbe(:,[2:end end]) - 4*phi_mbe;
f_lap2 = f_lap.^2;
f_x = (phi_mbe([1 1:end-1],:) - phi_mbe([2:end end],:))/2;
f_y = (phi_mbe(:,[1 1:end-1]) - phi_mbe(:,[2:end end]))/2;
f_g2 = f_x.^2 + f_y.^2;
f2 = (f_g2 - 1).^2;
energy2 = deta/2*sum(f_lap2(:)) + delta2/4*sum(f2(:));

disp(['1: ' num2str(energy) '   2: ' num2str(energy2)]);

end；